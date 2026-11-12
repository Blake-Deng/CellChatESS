# 十六进制追踪导出的整理说明

实际 C++ 直接追踪文件 ess_direct_trace.tsv 保留了正确的十六进制浮点字符串。
运行时 run_target_trace_CC22.R 第80行的 read.delim(trace_path) 自动将这些字符串转成数值；
因此旧 comparison_all_decisions.tsv 中 ESS_Pobs_hex / ESS_Pboot_hex 两列误导出成十进制。
这只是汇合表的导出类型问题，未影响实际追踪、严格比较、计数或 p 值。

本包不覆盖原运行文件。supplementary/S_G2_direct_decisions_cc16.tsv 和 cc22.tsv
用 reproduction/rebuild_trace_exports.py 从两侧原始直接追踪文件重新汇合，保留原始17位十进制和十六进制。
脚本验证每个浮点数十进制/十六进制可逆、delta、严格 > 判定及最终计数/p值，共500条/环境。
没有重新运行推断，也没有改变概率或 p 值公式。

未来执行入口的最小整理修复见 trace_export_colclasses.patch；该补丁尚未应用到本次实际运行脚本，
从而保留运行时源码身份。不能把重新汇合表说成一批新的实验。
