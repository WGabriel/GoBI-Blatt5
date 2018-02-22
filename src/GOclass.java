import java.util.HashMap;
import java.util.HashSet;

public class GOclass {
    public String id;
    public String name;
    public String namespace;
    // value = GOclass_id
    // is_a_updated is going to be "filled" with GO_id's from parents
    public HashSet<String> is_a_updated;
    public HashSet<String> is_a_original;
    // public HashSet<String> contains_GOs;
    // specifies distance to root (-1 = undefined) OPTIONAL
    public int level;
    public HashSet<String> associatedGenes;

    // Key = GO, String = Path (e.g. negative regulation of growth|response to stimulus|biological_process )
    public HashMap<String, String> shortestPathToGO;


    // [Term]
    // id: GO:0000001
    // name: mitochondrion inheritance
    // namespace: biological_process
    // def: "The distribution of mitochondria, including the mitochondrial genome,
    // into daughter cells after mitosis or meiosis, mediated by interactions
    // between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824,
    // PMID:11389764]
    // synonym: "mitochondrial inheritance" EXACT []
    // is_a: GO:0048308 ! organelle inheritance
    // is_a: GO:0048311 ! mitochondrion distribution

    public GOclass(String id, String name, String namespace) {
        this.id = id;
        this.name = name;
        this.namespace = namespace;
        this.is_a_updated = new HashSet<>();
        this.is_a_original = new HashSet<>();
        // this.contains_GOs = new HashSet<>();
        this.associatedGenes = new HashSet<>();
        this.level = -1;
        this.shortestPathToGO = new HashMap<>();
    }

}
