"""Prepare figure tables from the saved, version-verified CC22/ESS022 outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

B = Path(__file__).resolve().parent
D = B / 'figures/plot_data'
R = B / 'data/archived_plot_reference'
S = B / 'data/Figure1_DEJ_ESS022_existing_results_20261101'


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def write(table, name):
    table.to_csv(D / name, sep='\t', index=False, float_format='%.17g')


groups = read(S / 'group_counts.tsv')
legacy_groups = read(R / 'DE_groups.tsv')
assert groups.cells.is_unique and legacy_groups.cells.is_unique
mapping = groups.rename(columns={'group': 'cluster_id'}).merge(
    legacy_groups, on='cells', validate='one_to_one')
assert len(mapping) == 9 and int(mapping.cells.sum()) == 1000
rename = dict(zip(mapping.cluster_id, mapping.group))
write(mapping, 'DEJ_group_display_mapping.tsv')

weights = read(D / 'DE_network_weights.tsv')
weights['sender'] = weights.source_id.map(rename)
weights['receiver'] = weights.target_id.map(rename)
agreement = {}
for engine, name in [('Original CellChat', 'Figure1D_original_weights.tsv'),
                     ('CellChatESS', 'Figure1E_ess_weights.tsv')]:
    subset = weights[weights.engine.eq(engine)].copy()
    prior = read(R / name)
    joined = subset.merge(prior, on=['sender', 'receiver'], suffixes=('_new', '_old'),
                          validate='one_to_one')
    error = float(np.max(np.abs(joined.weight_new - joined.weight_old)))
    # Check the mapping against all 81 existing network weights, not counts alone.
    assert len(joined) == 81 and error < 1e-12
    agreement[engine + '_legacy_aggregate_max_abs_difference'] = error
    write(subset[['sender', 'receiver', 'weight', 'source_id', 'target_id']], name)

all_values = read(D / 'DEJ_nonzero_probabilities.tsv')
all_values['source'] = all_values.source_id.map(rename)
all_values['target'] = all_values.target_id.map(rename)
order = read(R / 'J_preserved_LR_order.tsv')
senders = read(R / 'J_preserved_sender_order.tsv')
assert len(order) == 16 and len(senders) == 3
new = all_values[all_values.interaction_name.isin(order.interaction_name)
                 & all_values.source.isin(senders.sender)
                 & all_values.probability.gt(0) & all_values.pvalue.lt(.05)].copy()
new = new.merge(order, on='interaction_name', validate='many_to_one')
prior_j = read(R / 'J_bubble_scores.tsv')
pairs = prior_j[['source', 'target', 'pair', 'pair_rank']].drop_duplicates()
present_pairs = new[['source', 'target']].drop_duplicates()
assert len(present_pairs.merge(pairs, on=['source', 'target'])) == len(present_pairs)
new = new.merge(pairs, on=['source', 'target'], validate='many_to_one')
lr = read(S / 'LRsig.tsv')
assert lr.interaction_name.is_unique
new = new.merge(lr.drop(columns=['lr_id']), on='interaction_name', validate='many_to_one')
new['visualization_score'] = -1 / np.log(new.probability)
assert np.isfinite(new.visualization_score).all()
# Keep compatibility columns with the original plot routine. Actual p-values
# remain in pvalue; the historical plotting category named pval is not reused.
new['prob'] = new.visualization_score
new['prob.original'] = new.probability
new = new.sort_values(['engine', 'lr_rank', 'pair_rank'])
write(new, 'J_bubble_scores.tsv')
for name in ['J_preserved_LR_order.tsv', 'J_preserved_sender_order.tsv', 'DE_groups.tsv']:
    write(read(R / name), name)

keys = ['engine', 'source', 'target', 'interaction_name']
overlap = new.merge(prior_j[keys + ['prob.original']], on=keys,
                    suffixes=('_new', '_legacy'), validate='one_to_one')
assert len(overlap) == len(new) == len(prior_j) == 146
agreement['legacy_bubble_max_abs_probability_difference'] = float(
    np.max(np.abs(overlap['prob.original_new'] - overlap['prob.original_legacy'])))
assert agreement['legacy_bubble_max_abs_probability_difference'] < 1e-14

original = new[new.engine.eq('Original CellChat')]
ess = new[new.engine.eq('CellChatESS')]
matched = original.merge(ess, on=['source', 'target', 'interaction_name'],
                         suffixes=('_original', '_ess'), validate='one_to_one')
assert len(matched) == 73 and matched.pvalue_original.equals(matched.pvalue_ess)
agreement.update({
    'network_edges_per_method': 81,
    'bubble_points_per_method': 73,
    'bubble_pairs': int(new.pair_rank.nunique()),
    'selected_LR': 16,
    'bubble_probability_max_abs_difference': float(np.max(np.abs(
        matched.probability_original - matched.probability_ess))),
    'bubble_pvalues_elementwise_exact': True,
    'inference_rerun': False,
    'display_mapping_basis': 'Nine distinct group sizes plus agreement with all 81 legacy network weights and all 146 legacy bubble records. Display annotations retained from the original figure; no new cell-type annotation.',
    'network_definition': 'Sum probabilities over all 1128 LR after setting p >= 0.05 to zero, matching aggregateNet default.',
    'bubble_definition': 'Original fixed set of 16 LR and three senders; probability > 0 and p < 0.05, re-evaluated on each new result.',
})
(B / 'DEJ_validation.json').write_text(json.dumps(agreement, indent=2) + '\n')
print(json.dumps(agreement, indent=2))
