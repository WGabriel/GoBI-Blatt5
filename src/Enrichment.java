public class Enrichment {
    // #GO:0099531
    // #GO:0040013
    // id fc signif
    // DNAJC25-GNG10 -1.3420 false
    // IGKV2-28 -2.3961 false
    // GSK3A -3.0193 true

    public String geneId;
    public Double fc;
    public boolean signif;
    // #: simulated enriched GO-ids
    // id: gene id
    // fc: log2 fold change estimation of gene diff expr (for KS-calculation)
    // signif=true, if gene was defined as diff expr. True is used to calculate
    // overrepres based enrich values
    // signif=false otherwise.

    public Enrichment(String geneId, Double fc, boolean signif) {
        this.geneId = geneId;
        this.fc = fc;
        this.signif = signif;
    }

}
