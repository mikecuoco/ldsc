#import "template.typ": preprint

#show: preprint.with(
  title: [Supplementary Information for\ ldsc-rs: Exact and approximate LD Score Regression at biobank scale],
  authors: (
    (name: "Haason, Sharif", affiliations: "1", corresponding: true),
    (name: "Khan, Yousef", affiliations: "2", corresponding: false),
  ),
  affils: (
    "1": "Independent",
    "2": "Independent",
  ),
  date: "March 2026",
  abstract: [
    This document collects supplementary tables referenced from the main
    text of ldsc-rs. Tables are numbered S1–S7 and each notes the main-text
    section it supports. The benchmark sweeps, complexity comparison, and
    per-subcommand feature inventory are placed here to keep the main-text
    table count focused on the headline results (parity, performance, the
    chunk-rounding finding, the cross-tool capability matrix, and the
    sketched-$hat(h)^2$ unbiasedness simulation).
  ],
  keywords: (
    [LD score regression],
    [heritability],
    [CountSketch],
    [supplementary],
  ),
  correspondence: "sharif@monark-markets.com",
)

// Supplementary tables number as S1, S2, … (main text uses 1, 2, …).
#set figure(numbering: n => "S" + str(n))

= Supplementary Tables

#figure(
  table(
    columns: 7,
    align: (left, right, right, right, right, right, right),
    table.header(
      [*$N$*],
      [*sketch-1000 + mask*\ *(h² truth cluster)*],
      [*sketch-1000 + mask*\ *peak RSS (MB)*],
      [*exact f32 wall (s)*],
      [*exact f32 RSS (MB)*],
      [*exact f64 wall (s)*],
      [*exact f64 RSS (MB)*],
    ),
    table.hline(stroke: 0.5pt),
    [503], [---#super("a")], [---], [9.4], [634], [15.5], [714],
    [10K], [7.4 s], [1,157], [93 s], [4,449], [243 s], [4,190],
    [20K], [6.4 s], [1,444], [133 s], [7,630], [482 s], [5,716],
    [50K], [9.7 s], [2,226], [292 s], [13,752], [---#super("b")], [---],
    [100K], [14.4 s], [4,022], [1{,}013 s], [14,609], [---#super("b")], [---],
  ),
  caption: [
    *Supports main-text §"Biobank-scale performance" (scaling figure).*
    Time and memory scaling across $N$ on a single Apple M5 Pro laptop
    (M5 Pro, 18 cores, 36 GB RAM, macOS 15). Wall time and peak RSS
    for the *h² truth-cluster sketch mode* (`--sketch 1000
    --snp-level-masking`), which recovers per-SNP exact heritability
    estimates within 0.001 across BMI and Height GWAS, against the
    two exact modes. At $N = 100{,}000$ the truth-cluster sketch
    delivers a $tilde$70#sym.times speedup over exact f32 (14.4 s
    vs.\ 1,013 s) and a $tilde$3.6#sym.times memory advantage (4.0 GB
    vs.\ 14.6 GB). At smaller $N$ the truth-cluster sketch is
    essentially indistinguishable in wall time from cheaper but
    biased sketch dimensions (e.g.\ `--sketch 200` finishes in
    13 s at $N = 100,000$ but gives h² off the truth cluster by
    ${tilde}0.013$ on real GWAS), so there is no speed reason to
    prefer a low $d$ over the recommended truth-cluster $d = 1000$.
    (#super("a")) `--sketch 1000` requires $d lt.eq N$ and is not
    applicable at $N = 503$ --- at 1000 Genomes scale, exact f32
    is the recommended mode. (#super("b")) Exact f64 at
    $N gt.eq 50,000$ takes multi-hour walls locally; the AWS EPYC
    7R13 figure of 727 s at $N = 50,000$ is reported in the
    main-text biobank performance table.
  ],
) <supp:scaling>

#figure(
  table(
    columns: 5,
    align: (left, right, right, right, right),
    table.header(
      [*$d$*], [*Measured $r$*], [*1000G time (s)*], [*Biobank time (s)*], [*Biobank speedup*],
    ),
    table.hline(stroke: 0.5pt),
    [50], [0.842], [3.6#super("†")], [16], [41#sym.times],
    [200], [0.960], [5.0#super("†")], [16], [41#sym.times],
    [500], [0.981], [7.9#super("†")], [18], [36#sym.times],
    [1000], [0.990], [10.5#super("†")], [20], [33#sym.times],
    [2000], [0.994], [28.3#super("†")], [22], [30#sym.times],
    [5000], [0.996], [---], [36], [18#sym.times],
    [10000], [0.997], [---], [59], [11#sym.times],
  ),
  caption: [
    *Supports main-text §"CountSketch accuracy–performance tradeoff".*
    CountSketch accuracy and runtime vs.\ sketch dimension. Accuracy ($r$):
    Pearson correlation with exact f64 LD scores, measured on the synthetic
    biobank dataset ($N = 50{,}000$, 1.66M SNPs; see main-text §Limitations).
    Times from the cost model $T(d) = T_"scatter" + T_"GEMM" dot d$ except where
    directly measured.
    #super("†")1000G local benchmarks (Ryzen 5 5600X).
    Biobank speedup vs.\ exact f64 (727 s).
  ],
) <supp:optimal-d>

#figure(
  table(
    columns: 5,
    align: (right, left, right, right, right),
    table.header(
      [*$N$*], [*Recommended $d$*], [*Wall (s)*],
      [*Pearson $r$*], [*$|hat(h)^2_("sketch") - hat(h)^2_("exact")|$*],
    ),
    table.hline(stroke: 0.5pt),
    [503], [500], [3.0], [0.994], [0.006],
    [10K], [1{,}000], [7.4], [0.994], [0.003],
    [20K], [1{,}000], [6.4], [0.993], [0.003],
    [50K], [1{,}000], [9.7], [0.994], [0.000],
    [100K], [1{,}000], [14.4], [0.993], [0.002],
  ),
  caption: [
    *Supports main-text §"Optimal $d$ as a function of $N$".*
    Sweet-spot $d$ per $N$ (criterion: smallest $d$ giving Pearson
    $r >= 0.993$ and $hat(h)^2$ within one regression-SE of exact),
    measured with `--snp-level-masking` on Apple M5 Pro. The optimum
    is approximately $N$-independent at $d approx 1{,}000$ because the
    sketch error is set by intra-window LD structure (window size
    $tilde$ 800-1500 SNPs), not by sample size. Wall scales mildly
    with $N$ from the fused $O(N c)$ scatter cost; at $N = 100{,}000$
    the recommended setting still finishes in $tilde$14 s, vs.\
    $tilde$17 min for exact f32 ($tilde$70 #sym.times speedup). For
    applications where exact per-SNP LD scores are critical
    (partitioned heritability, fine-mapping), bump to $d = 5{,}000$
    for $r >= 0.999$ at the cost of $tilde$30 s wall.
  ],
) <supp:optimal-d-by-N>

#figure(
  table(
    columns: 5,
    align: (right, right, right, right, right),
    table.header(
      [*`--chunk-size`*], [*mean L2*], [*Δ L2 vs c=50*], [*$hat(h)^2$ at h²=0.2*], [*$hat(h)^2$ at h²=0.5*],
    ),
    table.hline(stroke: 0.5pt),
    [25], [18.71], [$-0.13%$], [0.2115 (+5.8%)], [0.5580 (+11.6%)],
    [50 (Python default)], [18.73], [---], [0.2115 (+5.7%)], [0.5569 (+11.4%)],
    [100], [18.77], [+0.20%], [0.2109 (+5.4%)], [0.5556 (+11.1%)],
    [200 (ldsc-rs default)], [18.85], [+0.62%], [0.2107 (+5.3%)], [0.5528 (+10.6%)],
    [500], [18.91], [+0.95%], [0.2097 (+4.8%)], [0.5532 (+10.6%)],
    [1000], [19.05], [+1.68%], [0.2092 (+4.6%)], [0.5526 (+10.5%)],
  ),
  caption: [
    *Supports main-text §"The `--chunk-size` knob and implementation parity".*
    LD-score and heritability sensitivity to `--chunk-size`, computed on
    chromosome 22 of 1000 Genomes EUR ($N = 503$, 18,627 SNPs after MAF
    $gt.eq 0.05$, `--ld-wind-kb 1000`, `--global-pass`). Mean LD score grows
    monotonically with chunk size (range 1.8% of mean) due to chunked
    over-counting; $hat(h)^2$ correspondingly decreases. The full range of
    the chunk-size effect on $hat(h)^2$ is $tilde$1.1 percentage points of
    relative bias---small but systematic.
  ],
) <supp:chunksize-sweep>

#figure(
  table(
    columns: 4,
    align: (left, left, left, left),
    table.header(
      [*Aspect*], [*GCTA*], [*Python LDSC*], [*ldsc-rs*],
    ),
    table.hline(stroke: 0.5pt),
    [Window eviction], [Per-SNP exact], [c-SNP chunks, c=50], [c-SNP chunks, c=200 (default); per-SNP exact with `--snp-level-masking`],
    [Block boundaries], [Two-pass overlap averaging], [Single forward pass], [Single forward pass],
    [r² estimator default], [Biased ($+$noise floor)], [Unbiased], [Unbiased],
    [Unbiased r² flag], [`--ld-score-adj`], [n/a (always unbiased)], [n/a (always unbiased)],
    [Parallelism], [OpenMP per SNP], [Serial], [Rayon per chromosome + SIMD GEMM via faer],
    [Python-LDSC parity flag], [---], [---], [`--python-compat`],
  ),
  caption: [
    *Supports main-text §"Three implementations: GCTA, Python LDSC, ldsc-rs".*
    Algorithmic choices across the three independent LD-score implementations.
    All three share the same theoretical definition $ell_j = sum_k r^2_(j k)$
    but differ in windowing, edge handling, and r² estimator. GCTA's per-SNP
    exact semantics match the LDSC paper's mathematical statement; Python
    LDSC's chunked approximation is undocumented in the paper but inherited
    by every Python fork. ldsc-rs's `--snp-level-masking` reproduces GCTA's
    per-SNP semantics; `--python-compat` produces bit-identical Python LDSC
    output (`max_abs_diff = 0` verified on chr22 1000G EUR).
  ],
) <supp:implementations>

#figure(
  table(
    columns: 3,
    align: (left, center, left),
    table.header(
      [*Subcommand*], [*Status*], [*Notes*],
    ),
    table.hline(stroke: 0.5pt),
    [`munge-sumstats`], [#sym.checkmark], [Polars streaming; `--daner`/`--daner-n`; INFO-score handling (`--info`/`--info-min`) identical to Python LDSC],
    [`l2` (LD scores)], [#sym.checkmark], [`--sketch`, `--fast-f32`, `--snp-level-masking`],
    [`h2` (heritability)], [#sym.checkmark], [`--overlap-annot`, `--h2-cts`, two-step],
    [`rg` (genetic corr.)], [#sym.checkmark], [`--intercept-h2`, multi-trait],
    [`make-annot`], [#sym.checkmark], [BED interval annotation],
    [`cts-annot`], [#sym.checkmark], [Continuous annotation binning],
  ),
  caption: [
    *Supports main-text §"Feature parity".*
    Feature parity between ldsc-rs and Python LDSC. All subcommands and major
    flags are supported. ldsc-rs adds `--sketch`, `--fast-f32`, and
    `--snp-level-masking` as new modes not present in the original.
  ],
) <supp:features>

#figure(
  table(
    columns: 3,
    align: (left, left, left),
    table.header(
      [*Method*], [*Sketch cost*], [*GEMM cost*],
    ),
    table.hline(stroke: 0.5pt),
    [Exact], [---], [$O(N c w)$ per chunk],
    [Gaussian sketch], [$O(d N c)$], [$O(d c w)$],
    [CountSketch], [$O(N c)$], [$O(d c w)$],
  ),
  caption: [
    *Supports main-text Methods §"Complexity comparison".*
    Asymptotic cost comparison per chunk. $N$: individuals, $c$: chunk size,
    $w$: window size, $d$: sketch dimension. CountSketch achieves input
    sparsity time for the projection step. For PLINK BED genotypes (dense;
    $"nnz" approx N c$), the distinction between $O("nnz")$ and $O(N c)$
    vanishes, but the constant-factor advantage of scatter-add (one
    memory write per entry) over dense matrix--vector products ($d$
    multiply-accumulates per entry) remains critical. Gaussian sketch is
    shown for reference; ldsc-rs implements CountSketch only.
  ],
) <supp:sketch-complexity>
