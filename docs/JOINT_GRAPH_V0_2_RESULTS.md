# JG-0.2 frozen development results

The initial `JG-0.2` implementation and plan were frozen at commit `efdeee8`
before running 1,000 new-seed datasets. The full summary and replicate table
are retained in `research/joint_graph_v0_2_development_summary.tsv` and its
`_replicates.tsv` companion.

Most targeted repairs worked:

- matched null passed 50/50;
- off-grid moderate mediation passed 43/50;
- sparse pleiotropy passed 42/50;
- directional pleiotropy classified all 46 reportable replicates correctly;
- declared sample-overlap null passed 49/50;
- signed-LD null and mediation passed 47/50 and 44/50;
- correctly supplied fourfold scales passed 49/50;
- skeptical-prior mediation and diffuse-prior null passed 43/50 and 50/50.

The frozen run nevertheless failed its complete acceptance rule:

1. four directional and seven noisy-orientation analyses failed the original
   maximum-over-all-states evidence diagnostic;
2. mediation plus sparse pleiotropy passed only 29/50 with 20 independent
   blocks;
3. treating 70%-accurate orientation calls as certain produced high two-path
   support in 74.4% of reportable replicates.

The aligned boundary remained: exact aligned pleiotropy produced high two-path
support in 92% of replicates. Weak mediation produced high two-path support in
16%, demonstrating low and unstable power rather than a release-ready weak
signal regime.

`JG-0.2.1` addresses the first and third failures by exact integration over
orientation uncertainty and by restricting the evidence discrepancy gate to
posterior-relevant states. The coexistence result is treated as an information
limit: a new 30-block confirmatory cell is added while the failed 20-block cell
remains diagnostic. The `JG-0.2` seeds are development data and are not reused
for patch acceptance.
