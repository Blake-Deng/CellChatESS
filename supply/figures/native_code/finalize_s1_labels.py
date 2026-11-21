"""Apply the authors' final panel A and B labels to the archived Figure S1 PDF."""

import argparse
import hashlib
from itertools import groupby
from pathlib import Path

from pypdf import PdfWriter
from pypdf.generic import DecodedStreamObject


BASE_SHA256 = "69c8bf47223d3e295c2a9d189c4fe98c6b87f7bcb0c3fbe749b7b01fba307e22"


def label_line(text, x, y, widths, *, right=False):
    size = 8.55
    # The base PDF encodes an en dash at character zero in the second font subset.
    glyphs = [("/F2", "\0") if char == "–" else ("/F1", char) for char in text]
    if right:
        x -= sum(widths[font][ord(char)] * size / 1000 for font, char in glyphs)
    commands = ["q", f"1 0 0 1 {x:.5f} {y:.5f} cm", "BT", "0 0 Td"]
    for font, run in groupby(glyphs, key=lambda glyph: glyph[0]):
        encoded = "".join(char for _, char in run).encode("utf-16-be").hex()
        commands.extend([f"{font} {size} Tf", f"<{encoded}> Tj"])
    return "\n".join([*commands, "ET", "Q"]).encode("ascii")


def replace_two_lines(data, anchor, new_lines):
    start = data.rindex(b"q\n", 0, data.index(anchor))
    end_marker = b"\nET\nQ"
    first_end = data.index(end_marker, start) + len(end_marker)
    end = data.index(end_marker, first_end) + len(end_marker)
    return data[:start] + b"\n".join(new_lines) + data[end:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_pdf", type=Path)
    parser.add_argument("output_pdf", type=Path)
    args = parser.parse_args()
    if hashlib.sha256(args.base_pdf.read_bytes()).hexdigest() != BASE_SHA256:
        parser.error("Input must be the FigureS1.pdf archived at commit 67ddb4c.")
    if args.base_pdf.resolve() == args.output_pdf.resolve():
        parser.error("Choose a separate output path.")

    writer = PdfWriter(clone_from=args.base_pdf)
    page = writer.pages[0]
    widths = {}
    for key, font in page["/Resources"]["/Font"].items():
        values = iter(font.get_object()["/DescendantFonts"][0].get_object()["/W"])
        widths[key] = {}
        for start in values:
            for index, width in enumerate(next(values)):
                widths[key][start + index] = width
    # Keep the original rendering tool's newline after each source stream.
    data = b"".join(part.get_object().get_data() + b"\n" for part in page["/Contents"])
    data = replace_two_lines(data, b"240.4423973186 358.5666412157", [
        label_line("Median paired running-time ratio 27.48×", 428.60, 358.56664, widths, right=True),
        label_line("95% dataset-cluster CI 22.31–38.47×", 428.60, 348.30497, widths, right=True),
    ])
    data = replace_two_lines(data, b"509.153692275 539.9410349551", [
        label_line("Median of the 12 dataset-specific", 509.15369, 539.94103, widths),
        label_line("median paired running-time ratios 2.007×", 509.15369, 529.67536, widths),
        label_line("95% dataset-bootstrap CI 1.368–3.195×", 509.15369, 519.40969, widths),
    ])
    stream = DecodedStreamObject()
    stream.set_data(data)
    page.replace_contents(stream.flate_encode())
    writer.write(args.output_pdf)


if __name__ == "__main__":
    main()
