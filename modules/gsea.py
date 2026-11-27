import pandas as pd
from pydeseq2.ds import DeseqStats
import gseapy as gp

def create_gsea(dds, variant, control, gene_sets):
    s = DeseqStats(dds, contrast=["Type", variant, control])
    s.summary()
    res = s.results_df

    # remove empty stats, duplicates, then sort
    rnk = res.dropna(subset=["stat"])["stat"]
    rnk = rnk[~rnk.index.duplicated()]
    rnk = rnk.sort_values(ascending=False)

    if gene_sets == "MSigDB":
        gene_sets_name = "MSigDB_Hallmark_2020"
    elif gene_sets == "GO-BP":
        gene_sets_name = "GO_Biological_Process_2023"

    # compare genes to set of pathways
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets_name,
        min_size=5, max_size=500,
        permutation_num=1000,
        seed=42, threads=4,
        outdir=None,
    )

    # find the most significant up and down regulated pathways
    by_nes = pre.res2d.copy()
    by_nes["NES"] = by_nes["NES"].astype(float)
    by_nes = by_nes.sort_values("NES")
    plot_df = pd.concat([by_nes.head(10), by_nes.tail(10)]).drop_duplicates("Term")

    # plot the pathways
    for term in [by_nes.Term.iloc[0], by_nes.Term.iloc[-1]]:
        gp.gseaplot(rank_metric=pre.ranking, term=term, **pre.results[term])

    terms = by_nes.Term[:5].tolist() + by_nes.Term[-5:].tolist()
    gp.gseaplot2(
        terms=terms,
        RESs=[pre.results[t]["RES"] for t in terms],
        hits=[pre.results[t]["hits"] for t in terms],
        rank_metric=pre.ranking,
        figsize=(6, 8),
    )

    gp.dotplot(plot_df, column="FDR q-val", cutoff=1, top_term=len(plot_df),
               title=f"{variant} vs {control}", figsize=(4, 8),
               ofname=f"results/gsea_{variant}_{control}_{gene_sets}.png")
    