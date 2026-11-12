#!/usr/bin/env python3
"""Derive controlled research variants from the final candidate; no publication mutation."""
from pathlib import Path
import argparse,shutil,difflib,json,hashlib
P=argparse.ArgumentParser();P.add_argument('--output',type=Path,required=True);A=P.parse_args()
ROOT=Path(__file__).resolve().parents[1];SRC=ROOT/'software/CellChatESS'
if A.output.exists():P.error('Choose a fresh variants output directory')
A.output.mkdir(parents=True)
def replace_one(s,a,b):
 if s.count(a)!=1:raise ValueError('Non-unique/missing source anchor: '+a[:80])
 return s.replace(a,b,1)
base_cpp=(SRC/'src/cellchat_prob_cpp.cpp').read_text();base_r=(SRC/'R/accel_cellchat.R').read_text()
# A1 changes the projection only. Normalization still uses the original full matrix maximum.
start=base_r.index('  kernel_input <- subset_indexed_matrix(')
end=base_r.index('  data.use <- kernel_input$data',start)
a1=base_r[:start]+'''  kernel_input <- list(data = as.matrix(data.use), indices = list(
    ligand=ligand_idx, receptor=receptor_idx, coA=coA_idx, coI=coI_idx,
    agonist=agonist_idx, antagonist=antagonist_idx))
'''+base_r[end:]
# A4 evaluates every group and pair, including zero observed probabilities.
prefix,body=base_cpp.split('List cellchat_prob_ess_cpp(',1)
a4=replace_one(body,'if (L[k] > 0.0) active_senders.push_back(k);','active_senders.push_back(k);')
a4=replace_one(a4,'if (R[k] > 0.0) active_receivers.push_back(k);','active_receivers.push_back(k);')
a4=replace_one(a4,'if (val != 0.0) {','{ // Unpruned diagnostic comparator: retain every observed pair.')
a4=replace_one(a4,'pval[static_cast<R_xlen_t>(idx) + lr_offset] = static_cast<double>(reject[j]) / static_cast<double>(nboot);','pval[static_cast<R_xlen_t>(idx) + lr_offset] = pnull[idx] == 0.0 ? 1.0 : static_cast<double>(reject[j]) / static_cast<double>(nboot);')
# A3 further removes the cross-LR bootstrap cache, retaining per-LR scratch storage.
a3=replace_one(a4,'std::vector<double> tri_cache(static_cast<size_t>(U) * K * nboot, NA_REAL);','std::vector<double> tri_cache; // Cache disabled in A3.')
a3=replace_one(a3,'std::vector<unsigned char> tri_seen(static_cast<size_t>(U) * K * nboot, 0);','std::vector<unsigned char> tri_seen;')
lo=a3.index('        const int cache_gene = gene_to_cache[gene0];')
hi=a3.index('\n      }\n\n      for (int k : active_groups)',lo)
a3=a3[:lo]+'''        for (int k : active_groups) {
          bootAvg[local * K + k] = tri_mean_gene_cells(data, gene0, boot_cells[b * K + k]);
          tri_mean_evals_local += 1;
        }'''+a3[hi:]
a3=replace_one(a3,'Named("cache_slots") = static_cast<double>(static_cast<long long>(U) * static_cast<long long>(K) * static_cast<long long>(nboot)),','Named("cache_slots") = 0.0,')
# Resource-only derivative. Event format is compatible with run_xenium_profiled.R.
profile_helper=r'''
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <string>
static inline void emit_memory_profile_stage(const std::string& stage) {
 const char* path=std::getenv("CELLCHATESS_STAGE_FILE");
 if(path==NULL || path[0]=='\0')return;
 const auto now_ms=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
 std::ofstream out(path,std::ios::app);
 if(out)out << now_ms << '\t' << stage << "\t\t\t\t\t\t\n";
}
'''
profile=profile_helper+base_cpp
pp,pb=profile.split('List cellchat_prob_ess_cpp(',1)
pb=replace_one(pb,'  const int G = data.nrow();','  emit_memory_profile_stage("cpp_enter");\n  const int G = data.nrow();')
pb=replace_one(pb,'  NumericVector prob(Dimension(K, K, nLR));','  emit_memory_profile_stage("cpp_membership_cache_allocated");\n  NumericVector prob(Dimension(K, K, nLR));')
pb=replace_one(pb,'  return List::create(\n    Named("prob") = prob,','  emit_memory_profile_stage("cpp_probability_complete");\n  return List::create(\n    Named("prob") = prob,')
profile=pp+'List cellchat_prob_ess_cpp('+pb
# Trace-only derivative; input targets use zero-based LR and flattened pair indices.
trace_headers='\n#include <cstdlib>\n#include <fstream>\n#include <cstdio>\n#include <utility>\n'
tp,tb=base_cpp.split('List cellchat_prob_ess_cpp(',1)
trace_init=r'''
  const char* target_path=std::getenv("CELLCHATESS_TRACE_TARGETS");
  const char* output_path=std::getenv("CELLCHATESS_TRACE_OUTPUT");
  std::vector<std::pair<int,int> > trace_targets;
  if(target_path && target_path[0]) {
    std::ifstream input(target_path); if(!input)stop("Cannot open trace targets");
    int lr0,pair0; while(input>>lr0>>pair0)trace_targets.emplace_back(lr0,pair0);
  }
  std::ofstream trace_output;
  if(!trace_targets.empty()) {
    if(threads_used!=1)stop("Diagnostic tracing requires one thread");
    if(!output_path || !output_path[0])stop("Trace output path missing");
    trace_output.open(output_path,std::ios::app);if(!trace_output)stop("Cannot open trace output");
  }
  auto trace_decision=[&](int lr,int idx,int b,double observed,double boot) {
    if(std::find(trace_targets.begin(),trace_targets.end(),std::make_pair(lr,idx))==trace_targets.end())return;
    char line[1024];
    std::snprintf(line,sizeof(line),"%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%d\t%a\t%a\n",lr+1,idx,b+1,observed,boot,boot-observed,(boot-observed>0.0)?1:0,observed,boot);
    trace_output<<line;
  };
'''
tb=replace_one(tb,'  std::vector< std::vector<int> > boot_cells(nboot * K);',trace_init+'\n  std::vector< std::vector<int> > boot_cells(nboot * K);')
tb=replace_one(tb,'        if (val - pnull[idx] > 0.0) reject_local[j] += 1;','        trace_decision(lr,idx,b,pnull[idx],val);\n        if (val - pnull[idx] > 0.0) reject_local[j] += 1;')
trace=trace_headers+tp+'List cellchat_prob_ess_cpp('+tb
# Additional instrumentation locates output allocation/assignment boundaries. No model changes.
pp,pb=profile.split('List cellchat_prob_ess_cpp(',1)
pb=replace_one(pb,'  NumericVector pval(Dimension(K, K, nLR));',
 '  NumericVector pval(Dimension(K, K, nLR));\n  emit_memory_profile_stage("cpp_outputs_allocated");')
pb=replace_one(pb,'  for (int lr = 0; lr < nLR; ++lr) {\n    const bool simple_lr',
 '  emit_memory_profile_stage("cpp_lr_traversal_start");\n  for (int lr = 0; lr < nLR; ++lr) {\n    if(lr % 50 == 0) emit_memory_profile_stage(std::string("cpp_lr_") + std::to_string(lr));\n    const bool simple_lr')
profile=pp+'List cellchat_prob_ess_cpp('+pb
profile_r_helper = r"""
.ess_resource_stage <- function(stage) {
 path <- Sys.getenv('CELLCHATESS_STAGE_FILE','')
 if(!nzchar(path)) return(invisible(NULL))
 rss <- grep('^VmRSS:',readLines('/proc/self/status',warn=FALSE),value=TRUE)
 rss <- if(length(rss)) as.numeric(strsplit(trimws(sub('^VmRSS:','',rss[1])), '[[:space:]]+')[[1]][1]) else NA_real_
 cat(sprintf('%.0f\t%s\t%.0f\t\t\t\t\t\n',as.numeric(Sys.time())*1000,stage,rss),file=path,append=TRUE)
 invisible(NULL)
}
"""
profile_r=base_r
profile_r=replace_one(profile_r,'  res$n_genes_kernel <- n_genes_kernel','  .ess_resource_stage("R_after_cpp_return")\n  res$n_genes_kernel <- n_genes_kernel')
profile_r=replace_one(profile_r,'  dimnames(Pval) <- dimnames(Prob)','  dimnames(Pval) <- dimnames(Prob)\n  .ess_resource_stage("R_dimnames_assigned")')
profile_r=replace_one(profile_r,'  object@net <- list(prob = Prob, pval = Pval)','  object@net <- list(prob = Prob, pval = Pval)\n  .ess_resource_stage("R_outputs_assigned")')
profile_r=profile_r+profile_r_helper
variants={'A1_dense_full_gene' :(base_cpp,a1),'A3_no_cache_no_pruning':(prefix+'List cellchat_prob_ess_cpp('+a3,base_r),'A4_cache_no_pruning':(prefix+'List cellchat_prob_ess_cpp('+a4,base_r),'PROFILE':(profile,profile_r),'TRACE':(trace,base_r)}
summary=[]
for name,(cpp,r) in variants.items():
 dest=A.output/name/'CellChatESS';shutil.copytree(SRC,dest,ignore=shutil.ignore_patterns('*.so','*.o','*.dll','*.dylib','symbols.rds','.git','.DS_Store'))
 (dest/'src/cellchat_prob_cpp.cpp').write_text(cpp);(dest/'R/accel_cellchat.R').write_text(r)
 desc=(dest/'DESCRIPTION').read_text();(dest/'DESCRIPTION').write_text(desc.rstrip()+'\nConfig/CellChatESS/Build: '+name+'\n')
 patch=[]
 for rel in ['DESCRIPTION','R/accel_cellchat.R','src/cellchat_prob_cpp.cpp']:
  patch.extend(difflib.unified_diff((SRC/rel).read_text().splitlines(True),(dest/rel).read_text().splitlines(True),fromfile='final/'+rel,tofile=name+'/'+rel))
 (A.output/name/'controlled_diff.patch').write_text(''.join(patch))
 summary.append({'build':name,'release_version_base':'0.2.2','not_independent_release':True,'cpp_sha256':hashlib.sha256((dest/'src/cellchat_prob_cpp.cpp').read_bytes()).hexdigest(),'R_sha256':hashlib.sha256((dest/'R/accel_cellchat.R').read_bytes()).hexdigest()})
(A.output/'build_map.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Generated',len(summary),'controlled variants. A2 uses final algorithm=dense; A5 uses final algorithm=ess.')
