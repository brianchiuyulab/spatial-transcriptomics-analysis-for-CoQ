'''Export the existing Coq8a tissue maps with line-only scale bars.
No analysis, randomization, statistics, bead filters or coordinates are rerun.
The default source-map renderer continues to include its numerical scale label.
'''
import argparse, ast
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from render_proposal_spatial_maps import positive_barcodes
from spatial_plot_config import SAMPLES
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-dir',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    root=args.release_dir
    tree=ast.parse((root/'code/run_spatial_analysis.py').read_text(encoding='utf-8'))
    funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('plot_coq_map','add_scale_bar')]
    env=dict(globals());exec(compile(ast.Module(body=funcs,type_ignores=[]),'source_renderers','exec'),env)
    plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.labelsize':12,'axes.titlesize':15,'axes.linewidth':1})
    meta=pd.read_csv(root/'data/GSE266933_compact_metadata.csv.gz')
    for sample in ('Young','Geriatric'):
        frame=meta[meta['sample']==sample].copy()
        _,positive=positive_barcodes(root/'data'/SAMPLES[sample]['file'])
        frame['Coq8a_detected']=frame.barcode.isin(positive)
        env['plot_coq_map'](frame,sample,args.output/f'{sample}_baronly.png',scale_label=False)
if __name__=='__main__':main()
