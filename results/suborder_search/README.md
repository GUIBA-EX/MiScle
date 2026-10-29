# M suborder search / M亚目子引物搜索

All M endpoints are fixed at F1561–1580/R1741–1760 on NC_003522.1 ungapped 16S. Ordinary and D mixtures each contain one shared F and two R oligos. / M所有版本参考端点相同；普通与D各使用一条共用F及两条R。

- `strict_selected_mixtures.json`: original ≤20 strict WC/GT and ≤100 WC3/GT2 screen; it did not yield a useful Refertina subprimer. / 原低背景筛选缺少充分Refertina补偿。
- `low_background_selected_mixtures.json`: ≤50 WC, ≤150 GT and ≤800 WC3; ordinary union 282 covered/13 exact, D union 284/19, WC/GT 27/111. / 较低WC方案完全匹配较低。
- `selected_mixtures.json`: ≤150 WC/GT, ≤1200 WC3; ordinary union 282/78, D union 284/98, WC/GT 103/111. The increased WC count buys exact matches, not additional species coverage. / 当前选择增加WC背景以换取完全匹配，不增加相对较低WC方案的覆盖。
- Within two species of the best union coverage, maximize exact matches, then coverage and the weaker suborder's coverage. This is an exploratory heuristic search. / 在最优覆盖损失≤2物种范围内优先完全匹配，为探索性启发式筛选。

Install `primer3-py` for candidate generation and thermal selection. The following commands run from the repository root, using a new empty working directory (for example `/tmp/miscle-reproduction`):

```sh
python3 scripts/generate_suborder_candidates.py . /tmp/miscle-reproduction
c++ -O3 -std=c++17 scripts/scan_suborder_sites.cpp -o /tmp/miscle-sites
/tmp/miscle-sites /tmp/miscle-reproduction/target/sites
/tmp/miscle-sites /tmp/miscle-reproduction/background/sites
python3 scripts/choose_suborder_mixtures.py /tmp/miscle-reproduction .
```

Frozen candidate pools, selected mixtures, raw pair hits, independent verification inputs and logs are retained here. / 本目录保留候选池、组合选择、原始命中、独立复核输入及日志。

`python3 scripts/verify.py` independently rechecks the frozen final and cross-mixture paired results with the Python standard library. / 独立复核无需安装额外Python库。

Rebuild release tables and the six base-matching SVGs after verification: `python3 scripts/build_suborder_release.py`. The reference-position plot uses `python3 scripts/plot_regions.py` with matplotlib installed. / 复核后可重建统计表与六张碱基匹配图，参考位置图另用matplotlib绘制。
