# Independent suborder F/R search / 亚目独立F/R搜索

M endpoints are fixed at F1561–1580/R1741–1760 on ungapped NC_003522.1 16S. Refertina and Vacatina each have independent F and R candidate pools; ordinary and D endpoints coincide. / 固定M坐标，分别筛选两个亚目的F/R，普通和D端点一致。

The search contains 105 unique oligos. Each D oligo has degeneracy ≤8. Maximum possible homopolymer length is five, matching the Primer3 default; concrete variants with six consecutive identical bases are excluded. / 105条独立寡核苷酸，每条简并≤8，允许最长五连、排除可能六连。

- `selected.json`: released suborder pairs / 正式子引物。
- `selected_M_only.json`: independent M-only optimum before the lower-background R-D choice / 降低Refertina背景前的M单窗口对照。
- `selected_delta5.json`: ≤5°C midpoint-difference comparator; ordinary MV exactly matches one Vacatina / 普通MV严格Tm差分支。
- `screened_pairs.json`: paired target/background counts and thermal evidence / 配对覆盖、背景与热力学证据。
- `verification_selected.json`: six final plus sixteen optional same-tube extra pairs / 六个正式配对与十六个潜在同管额外配对。

Per-suborder candidates within two species of the best coverage are ranked by exact matching, coverage and background. The released R-D option retains the best SD+M coverage and exact union, then minimizes the higher background model count, considers Refertina exact matching and minimizes total background. / 亚目内在最佳覆盖损失≤2物种范围内优先完全匹配；正式R-D在SD+M覆盖及完全匹配并集不降低的条件下压低背景，并继续考虑Refertina完全匹配及背景总数。

The same-tube four-oligo model is a sensitivity check, not the recommended PCR layout. Released combinations describe separate suborder PCRs. / 四条引物同管只作敏感性检查，正式组合按亚目分管PCR解释。

From the repository root, install `requirements-design.txt` to regenerate candidates and use a new empty directory:

```sh
python3 scripts/generate_independent_candidates.py . /tmp/miscle-independent-reproduction
c++ -O3 -std=c++17 scripts/scan_suborder_sites.cpp -o /tmp/miscle-sites
/tmp/miscle-sites /tmp/miscle-independent-reproduction/target/sites
/tmp/miscle-sites /tmp/miscle-independent-reproduction/background/sites
python3 scripts/choose_independent_pairs.py /tmp/miscle-independent-reproduction .
python3 scripts/prepare_independent_verification.py . /tmp/miscle-independent-reproduction
c++ -O3 -std=c++17 scripts/independent_scan.cpp -o /tmp/miscle-pairs
/tmp/miscle-pairs /tmp/miscle-independent-reproduction/verify
cp /tmp/miscle-independent-reproduction/verify/independent_hits.tsv /tmp/miscle-independent-reproduction/independent_hits.tsv
python3 scripts/verify.py . /tmp/miscle-independent-reproduction
```

To verify the frozen release without installing design dependencies or using a network, run `python3 scripts/verify.py`. / 仅复核冻结发布结果时无需安装设计依赖或联网。

`python3 scripts/build_suborder_release.py` rebuilds release tables and six base-matching SVGs from these frozen scans. / 从冻结扫描重建统计和六张碱基匹配图。
