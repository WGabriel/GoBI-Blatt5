import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.*;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

public class Runner {

    public static void main(String[] args) {
        File obo = new File(""); // GO obo_file
        String root = ""; // GO namespace [biological_process, cellular_component, molecular_function]
        File mapping = new File(""); // gene2go_mapping
        String mappingtype = ""; // [ensembl|go] --> defines how mapping is processed
        boolean overlapoutGiven = false;
        File overlapout = new File(""); // NOT OBLIGATORY overlap_out_tsv
        File enrich = new File(""); // simulation file
        File output = new File(""); // output_tsv
        int minsize = 0;
        int maxsize = 0;

        if (args.length < 16 || args.length > 18) {
            System.err.println("Required 16-18 parameters. Found: " + args.length);
            StringBuilder consoleParameters = new StringBuilder();
            for (int i = 0; i < args.length; i++) {
                consoleParameters.append(" args[").append(i).append("]=").append(args[i]).append("\n");
            }
            System.err.println("ConsoleParameters:" + consoleParameters);
        } else {
            // Console parameters are correct!
            for (int i = 0; i < args.length; i = i + 2) {
                switch (args[i]) {
                    case "-obo":
                        obo = new File(args[i + 1]);
                        continue;
                    case "-root":
                        root = args[i + 1];
                        if (!root.equals("biological_process") && !root.equals("cellular_component")
                                && !root.equals("molecular_function")) {
                            System.err.println(
                                    "The parameter root does not equal biological_process/cellular_component/molecular_function!");
                        }
                        continue;
                    case "-mapping":
                        mapping = new File(args[i + 1]);
                        continue;
                    case "-mappingtype":
                        mappingtype = args[i + 1];
                        if (!mappingtype.equals("go") && !mappingtype.equals("ensembl")) {
                            System.err.println("The parameter mappingtype does not equal go/ensembl!");
                        }
                        continue;
                    case "-overlapout":
                        overlapoutGiven = true;
                        overlapout = new File(args[i + 1]);
                        continue;
                    case "-enrich":
                        enrich = new File(args[i + 1]);
                        continue;
                    case "-o":
                        output = new File(args[i + 1]);
                        continue;
                    case "-minsize":
                        minsize = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-maxsize":
                        maxsize = Integer.parseInt(args[i + 1]);
                }
            }
            System.out.print("|UserInput|");
            System.out.println("obo:\t" + obo.getAbsolutePath());
            System.out.println("root:\t" + root);
            System.out.println("mapping:\t" + mapping.getAbsolutePath());
            System.out.println("mappingtype:\t" + mappingtype);
            if (overlapoutGiven) {
                System.out.println("overlapout:\t" + overlapout.getAbsolutePath());
            } else {
                System.out.println("overlapout:\tNOT GIVEN");
            }
            System.out.println("enrich:\t" + enrich.getAbsolutePath());
            System.out.println("output:\t" + output.getAbsolutePath());
            System.out.println("minsize:\t" + minsize);
            System.out.println("maxsize:\t" + maxsize);
        }

        // namespaces = library with namespaces and all GO-terms parset out of obo-file
        HashMap<String, HashMap<String, GOclass>> namespaces = parseOboFile(obo);

        // fill is_a_updated with GO_id's from parents
        namespaces = updateGOlinks(namespaces, root);

//        String s = getStringShortestStringPath("GO:0006641", "GO:0044238", namespaces.get(root));
//        System.out.println("PATH: GO:0006641->GO:0044238: " + s);
//        System.out.println("INVERTED: " + invertPath(s));

        // parse associated gene_ids of GOs from mapping file
        namespaces = parseAssociatedGeneIds(mappingtype, mapping, namespaces, root);
        // The hashSet associatedGenes has to be updated (= add genes of children)
        namespaces = updateAssociatedGeneIds(namespaces, root);

//        File completeNamespace = new File(
//                "C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt5\\GOEnrichment\\result\\completeRootNamespace.txt");
//        writeNamespaces(namespaces, root, completeNamespace);

        if (overlapoutGiven) {
            createOverlapout(namespaces, root, overlapout, minsize, maxsize);
        }

        // parse simul_exp_go_bp_ensembl.tsv
        HashMap<String, Enrichment> enrichGenes = parseEnrichFileGenes(enrich);
        HashSet<String> enrichGOs = parseEnrichFileGO(enrich);

        // creates tsv file containing the resulting enrichment info
        createEnrichmentOutput(namespaces.get(root), minsize, maxsize, enrichGenes, enrichGOs, output);

        System.out.println("End of main.");
    }

    private static void createEnrichmentOutput(HashMap<String, GOclass> namespace, int minsize, int maxsize,
                                               HashMap<String, Enrichment> enrichGenes, HashSet<String> enrichGOs, File output) {
        System.out.println("Begin: createEnrichmentOutput...");
        ArrayList<Object[]> result = new ArrayList<>();
        Object[] header = new Object[1];
        header[0] = "term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true";
        result.add(header);

        // pop_size = amount of genes in root AND in enrichment file
        // --> intersection. Required in HG-calculation.
        HashSet<String> genesRootIntersectEnrich = new HashSet<>();
        for (GOclass go : namespace.values()) {
            genesRootIntersectEnrich.addAll(go.associatedGenes);
        }
        genesRootIntersectEnrich.retainAll(enrichGenes.keySet());

        int population_size = genesRootIntersectEnrich.size();

        // nrOfSuccesses = Amount of signif genes in popSize
        int nrOfSuccesses = 0;
        for (String gene : genesRootIntersectEnrich) {
            if (enrichGenes.containsKey(gene)) {
                if (enrichGenes.get(gene).signif) {
                    nrOfSuccesses++;
                }
            }
        }

        for (GOclass go : namespace.values()) {
            if (go.associatedGenes.size() >= minsize && go.associatedGenes.size() <= maxsize) {
                String term = go.id;
                String name = go.name;
                //System.out.println("*****Creating Result Line for term: " + term + " name: " + name);
                int size = 0;
                boolean is_true = false;
                int noverlap = 0;
                double hg_pval;
                double hg_fdr = 0;
                double fej_pval;
                double fej_fdr = 0;
                double tks_stat;
                double tks_pval;
                double tks_fdr = 0;
                String shortest_path_to_a_true = "";

                // calculate 'size', 'noverlap' and 'inSetFoldChangeDistribution' (for KS)
                // key = geneId, value = fold change
                HashMap<String, Double> inSetFoldChangeDistr = new HashMap<>();
                for (String assGene : go.associatedGenes) {
                    if (enrichGenes.containsKey(assGene)) {
                        size++;
                        inSetFoldChangeDistr.put(assGene, enrichGenes.get(assGene).fc);
                        if (enrichGenes.get(assGene).signif) {
                            noverlap++;
                        }
                    }
                }

                // calculate 'is_true'
                if (enrichGOs.contains(go.id)) {
                    is_true = true;
                }

                // ***hg_vpal***
                // population size = # genes in enrichtment-file + root
                // nr of successes = # genes in population, where signinf=true
                // sample size = 'size' (from assigment)
                HypergeometricDistribution hg = new HypergeometricDistribution(population_size, nrOfSuccesses, size);
                hg_pval = hg.upperCumulativeProbability(noverlap);
                // ***hg_fdr ***--> calculated below

                // ***fej_pval***
                hg = new HypergeometricDistribution(population_size - 1, nrOfSuccesses - 1, size - 1);
                fej_pval = hg.upperCumulativeProbability(noverlap - 1);

                //***KolmogorovSmirnov***
                // create a copy of enrichGenes containing only complementary set of enrichGenes and inSetFoldChangeDistr
                HashMap<String, Enrichment> backgroundFoldChangeDistr = new HashMap<>(enrichGenes);
                // backgroundFoldChangeDistr complementary retaining:
                for (Iterator<Map.Entry<String, Enrichment>> it = backgroundFoldChangeDistr.entrySet().iterator();
                     it.hasNext(); ) {
                    Map.Entry<String, Enrichment> entry = it.next();
                    if (!genesRootIntersectEnrich.contains(entry.getKey())) {
                        it.remove();
                    } else if (inSetFoldChangeDistr.containsKey(entry.getKey())) {
                        it.remove();
                    }
                }

                // convert inSetFoldChangeDistr and backgroundFoldChangeDistr to arrays
                double[] inSetFCArr = new double[inSetFoldChangeDistr.size()];
                int counter = 0;
                for (double d : inSetFoldChangeDistr.values()) {
                    inSetFCArr[counter] = d;
                    counter++;
                }

                //copy fc values of backgroundFoldChangeDistr to array
                double[] bgFCArr = new double[backgroundFoldChangeDistr.size()];
                counter = 0;
                for (Enrichment e : backgroundFoldChangeDistr.values()) {
                    bgFCArr[counter] = e.fc;
                    counter++;
                }

                KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
                // ks_stat in_set_distribution , background-distribution-array = alle ohne in_set
                tks_stat = ks.kolmogorovSmirnovStatistic(inSetFCArr, bgFCArr);
                tks_pval = ks.kolmogorovSmirnovTest(inSetFCArr, bgFCArr);

                if (is_true || enrichGOs.isEmpty()) {
                    //assignment: "shortest_path_to_a_true is empty, if is_true=true or no true entries provided"
                } else {
                    shortest_path_to_a_true = getShortestPathToTrue(term, enrichGOs, namespace);
                }

                Object[] line = new Object[13];
                line[0] = term;
                line[1] = name;
                line[2] = size;
                line[3] = is_true;
                line[4] = noverlap;
                line[5] = hg_pval;
                line[6] = hg_fdr;
                line[7] = fej_pval;
                line[8] = fej_fdr;
                line[9] = tks_stat;
                line[10] = tks_pval;
                line[11] = tks_fdr;
                line[12] = shortest_path_to_a_true;
                result.add(line);
                // result.add(term + "\t" + name + "\t" + size + "\t" + is_true + "\t" + noverlap + "\t" + hg_pval + "\t" + hg_fdr + "\t" + fej_pval + "\t" + fej_fdr + "\t" + tks_stat + "\t" + tks_pval + "\t" + tks_fdr + "\t" + shortest_path_to_a_true);
            }
        }

        //needed for benjamini-hochberg-correction
        Double[] hg_pvalues = new Double[result.size() - 1];
        Double[] fej_pvalues = new Double[result.size() - 1];
        Double[] ks_pvalues = new Double[result.size() - 1];
        for (int i = 1; i < result.size(); i++) {
            hg_pvalues[i - 1] = (Double) result.get(i)[5];
            fej_pvalues[i - 1] = (Double) result.get(i)[7];
            ks_pvalues[i - 1] = (Double) result.get(i)[10];
        }

        double[] hg_fdr_values = getBenjaminiHochbergCorrection(hg_pvalues);
        double[] fej_fdr_values = getBenjaminiHochbergCorrection(fej_pvalues);
        double[] ks_fdr_values = getBenjaminiHochbergCorrection(ks_pvalues);

        // convert Object-Array (result) to String-Array (resultString) for PrintOut
        ArrayList<String> resultString = new ArrayList<>();
        // add header to resultString
        resultString.add((String) result.get(0)[0]);
        for (int i = 1; i < result.size(); i++) {
            StringBuilder tempLine = new StringBuilder();
            // update hg_fdr with adjusted pValues
            result.get(i)[6] = hg_fdr_values[i - 1];
            // update fej_fdr with adjusted pValues
            result.get(i)[8] = fej_fdr_values[i - 1];
            // update ks_fdr with adjusted pValues
            result.get(i)[11] = ks_fdr_values[i - 1];
            for (Object o : result.get(i)) {
                if (o instanceof Double) {
                    o = String.format("%g", o);
                }
                tempLine.append(o).append("\t");
            }
            resultString.add(tempLine.toString());
        }
        // write resultString to output-location
        writeArrayList(resultString, output);
        System.out.println("Finished: createEnrichmentOutput.");
    }

    public static String getShortestPathToTrue(String goId, HashSet<String> enrichGOs, HashMap<String, GOclass> namespace) {
        // System.out.println("|getShortestPathToTrue| Begin calculation with " + goId);

        // Add paths GO->CommonGO plus enrichGO->CommonGO
        ArrayList<String> shortestPathGoToCommon = new ArrayList<>();
        ArrayList<String> shortestEnrichGOToCommon = new ArrayList<>();
        double shortestPathLength = Integer.MAX_VALUE;

        for (String goId_enrich : enrichGOs) {
            // For every geneId in enrichGOs, gather common GO's
            // System.out.println("Current enrichGO: " + goId_enrich);

            //error checking
            if (goId_enrich.equals(goId)) {
                System.err.println("|getShortestPathToTrue| goId is not allowed to occur in enrichGO's: " + goId);
            }
            //commonGOs is an INTERSECTION of enrich_goId.is_a and goId_is_a
            HashSet<String> commonGOs = new HashSet<>(namespace.get(goId).is_a_updated); // use the copy constructor
            commonGOs.retainAll(namespace.get(goId_enrich).is_a_updated);

            // add goId, if it occurs in goId_enrich.is_a_updated
            if (namespace.get(goId_enrich).is_a_updated.contains(goId)) {
                commonGOs.add(goId);
            }
            // add goId_enrich, if it occurs in goId.is_a_updated
            if (namespace.get(goId).is_a_updated.contains(goId_enrich)) {
                commonGOs.add(goId_enrich);
            }
            // End: gathering common GO's

            if (!commonGOs.isEmpty()) {
                for (String common : commonGOs) {
                    // constructing path GO-->Common
                    ArrayList<String> tempPathGoToCommon = new ArrayList<>();
                    if (goId.equals(common)) {
                        System.out.println("goId: " + goId + " common: " + common);
                        tempPathGoToCommon = new ArrayList<>();
                        tempPathGoToCommon.add(namespace.get(goId).name);
                    } else {
                        try {
                            tempPathGoToCommon = getShortestStringPath(goId, common, namespace);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }

                    // constructing path EnrichGO-->Common
                    ArrayList<String> tempEnrichGOToCommon = new ArrayList<>();
                    if (goId_enrich.equals(common)) {
                        System.out.println("enrich_goId: " + goId_enrich + " common: " + common);
                        tempEnrichGOToCommon = new ArrayList<>();
                        tempEnrichGOToCommon.add(namespace.get(goId_enrich).name);
                    } else {
                        try {
                            tempEnrichGOToCommon = getShortestStringPath(goId_enrich, common, namespace);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                    if ((tempPathGoToCommon.size() + tempEnrichGOToCommon.size()) < shortestPathLength) {
                        shortestPathLength = (tempPathGoToCommon.size() + tempEnrichGOToCommon.size());

                        shortestPathGoToCommon.clear();
                        shortestPathGoToCommon = (ArrayList<String>) tempPathGoToCommon.clone();
                        shortestEnrichGOToCommon.clear();
                        shortestEnrichGOToCommon = (ArrayList<String>) tempEnrichGOToCommon.clone();
                    }
                }
            } else {
                // if there is no common GO --> use root (?)
                System.err.println("|getShortestPathToTrue| No common GO with " + goId);
            }
        }

        //convert ArrayList<String> shortestPathGoToCommon to result String
        String result = "";
        for (int i = 0; i < shortestPathGoToCommon.size(); i++) {
            if (i == 0) {
                result = shortestPathGoToCommon.get(i);
            } else {
                result = result.concat("|").concat(shortestPathGoToCommon.get(i));
            }
        }
        result = result.concat(" * ");
        //convert ArrayList<String> shortestEnrichGOToCommon to result String
        if (!shortestEnrichGOToCommon.get(shortestEnrichGOToCommon.size() - 1).equals(shortestPathGoToCommon.get(shortestPathGoToCommon.size() - 1))) {
            System.err.println("Go_id: " + goId + " shortestPathGoToCommon: " + shortestPathGoToCommon.toString() + " and shortestEnrichGOToCommon " + shortestEnrichGOToCommon.toString());
            System.err.println("letztes shortestPathGoToCommon: " + shortestPathGoToCommon.get(shortestPathGoToCommon.size() - 1) + " != letztes shortestEnrichGOToCommon: " + shortestEnrichGOToCommon.get(shortestEnrichGOToCommon.size() - 1));
            // result = result.concat("|").concat(shortestEnrichGOToCommon.get(shortestEnrichGOToCommon.size() - 1));
        }

        for (int i = shortestEnrichGOToCommon.size() - 2; i >= 0; i--) {
            // invert Path of tempEnrichGOToCommon
            // Note: tempEnrichGOToCommon.size() - 2 to avoid duplicate
            result = result.concat("|").concat(shortestEnrichGOToCommon.get(i));
        }

        // System.out.println("|getShortestPathToTrue| Shortest Path for " + goId + " to enrichGOs is: " + result);
        return result;
    }

    public static ArrayList<String> getShortestStringPath(String childGoId, String parentGoId, HashMap<String, GOclass> namespace) throws Exception {
        if (childGoId.equals(parentGoId)) {
            System.err.println("|getShortestStringPath| childGoId must be different from parentGoId. (" + childGoId + ")");
            throw new Exception("This path leads not to a result.");
        } else if (namespace.get(childGoId).shortestPathToGO.containsKey(parentGoId)) {
            // Avoid double calculation
            return namespace.get(childGoId).shortestPathToGO.get(parentGoId);
        } else if (namespace.get(childGoId).is_a_original.contains(parentGoId)) {
            ArrayList<String> result = new ArrayList<>();
            result.add(namespace.get(childGoId).name);
            result.add(namespace.get(parentGoId).name);
            namespace.get(childGoId).shortestPathToGO.put(parentGoId, result);
            if ((!result.get(0).equals(namespace.get(childGoId).name)) || (!result.get(result.size() - 1).equals(namespace.get(parentGoId).name))) {
                // Error checking
                System.err.println("|getShortestStringPath| is_a_original.contains(parentGo) First/last result does not equal child/parent. First: " + result.get(0) + " (" + namespace.get(childGoId).name + ") Last: " + result.get(result.size() - 1) + " (" + namespace.get(parentGoId).name + ")");
            }
            return result;
        } else if (namespace.get(childGoId).is_a_original.isEmpty()) {
            // This path leads not to a result
            throw new Exception("This path leads not to a result.");
        } else {
            // Perform 'usual' calculation
            double bestLength = Integer.MAX_VALUE;
            ArrayList<String> result = new ArrayList<>();

            for (String go : namespace.get(childGoId).is_a_original) {
                ArrayList<String> tempResult = new ArrayList<>();
                try {
                    tempResult = getShortestStringPath(go, parentGoId, namespace);
                } catch (Exception e) {
                    // e.printStackTrace();
                }

                if (tempResult.size() < bestLength /*&& tempResult.size() > 0*/) {
                    result.clear();
                    result.add(namespace.get(childGoId).name);
                    result.addAll((ArrayList<String>) tempResult.clone());
                    bestLength = tempResult.size();
                }
            }
            if (result.size() < 1) {
                System.err.println("|getShortestStringPath| Best result cannot be < 1. child: " + childGoId + " parent: " + parentGoId);
            }
            if (!result.get(result.size() - 1).equals(namespace.get(parentGoId).name)) {
                // Add parentGoId.name to the last recursive iteration
                result.add(namespace.get(parentGoId).name);
            }
            namespace.get(childGoId).shortestPathToGO.put(parentGoId, result);
            if ((!result.get(0).equals(namespace.get(childGoId).name)) || (!result.get(result.size() - 1).equals(namespace.get(parentGoId).name))) {
                // Error checking
                System.err.println("|getShortestStringPath| 'usual' First/last result does not equal child/parent. First: " + result.get(0) + " (" + namespace.get(childGoId).name + ") Last: " + result.get(result.size() - 1) + " (" + namespace.get(parentGoId).name + ")");
            }
            return result;
        }
    }


    public static String invertPath(String path) {
        //Invert GO:0001|GO:0002|GO:0003 to GO:0003|GO:0002|GO:0001
        String[] originalOrder = path.split("\\|");
        String result = "";
        for (int i = originalOrder.length - 1; i >= 0; i--) {
            result = result.concat(originalOrder[i] + "|");
        }
        return result.substring(0, result.length() - 1);
    }

    public static double[] getBenjaminiHochbergCorrection(Double[] pValuesAsArray) {
        //Convert pValuesAsArrayToList
        List<Double> pValues = Arrays.asList(pValuesAsArray);

        List<Double> sortedPValues = new LinkedList<>(pValues);
        Collections.sort(sortedPValues); // ascending

        int n = pValues.size();
        double[] sortedAdjustedPvalues = new double[n];
        double[] adjustedPvalues = new double[n];
        double mFDR = sortedPValues.get(n - 1);
        sortedAdjustedPvalues[n - 1] = mFDR;

        for (int k = n - 1; k > 0; k--) {
            mFDR = Math.min(sortedPValues.get(k - 1) * ((double) n / k), mFDR);
            sortedAdjustedPvalues[k - 1] = mFDR;
        }

        for (int i = 0; i < pValues.size(); i++) {
            adjustedPvalues[i] = sortedAdjustedPvalues[sortedPValues.indexOf(pValues.get(i))];
        }

        return adjustedPvalues;
    }

    public static HashSet<String> parseEnrichFileGO(File enrich) {
        System.out.println("Begin: parseEnrichFileGO...");
        // list of GOs in enrich file
        HashSet<String> result = new HashSet<String>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(enrich));
            String line = null;
            while ((line = br.readLine()) != null) {
                // #GO:0099531
                if (line.startsWith("#")) {
                    result.add(line.substring(1, line.length()));
                    // System.out.println("#: " + line.substring(1, line.length()));
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished: parseEnrichFileGO. GOs in Enrich-Output: " + result.size());
        return result;
    }

    public static HashMap<String, Enrichment> parseEnrichFileGenes(File enrich) {
        System.out.println("Begin: parseEnrichFileGenes...");
        // key = geneId, value = Enrichment
        HashMap<String, Enrichment> result = new HashMap<String, Enrichment>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(enrich));
            String line = null;
            while ((line = br.readLine()) != null) {
                // #GO:0099531
                // id fc signif
                // DNAJC25-GNG10 -1.3420 false
                // GSK3A -3.0193 true
                if (line.startsWith("#")) {
                    // see method parseEnrichFileGOs
                } else if (line.equals("id\tfc\tsignif")) {
                    // ignore header (id fc signif)
                    // System.out.println("HEADER: " + line);
                } else {
                    String geneId = line.split("\t")[0];
                    Double fc = Double.parseDouble(line.split("\t")[1]);
                    boolean signif = Boolean.parseBoolean(line.split("\t")[2]);
                    if (geneId.isEmpty() || line.split("\t").length < 3) {
                        System.err.println("Error in line: " + line);
                    }
                    if (result.containsKey(geneId)) {
                        System.err.println("geneId: " + geneId + " already exists in enrichment HashMap.");
                    } else {
                        result.put(geneId, new Enrichment(geneId, fc, signif));
                        // System.out.println("geneId: " + geneId + " fc: " + fc + " signif: " +
                        // signif);
                    }
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished: parseEnrichFileGenes. GeneIds in Enrich-Output: " + result.size());
        return result;
    }

    private static HashMap<String, HashMap<String, GOclass>> updateAssociatedGeneIds(HashMap<String, HashMap<String, GOclass>> namespaces, String root) {
        System.out.println("Begin: updateAssociatedGeneIds...");
        // Add the associated genes of the children to every GOclass
        for (GOclass go : namespaces.get(root).values()) {
            // iterate through all GOs in namespace
            for (String GOclass_id : go.is_a_updated) {
                namespaces.get(root).get(GOclass_id).associatedGenes.addAll(go.associatedGenes);
            }
        }
        System.out.println("Finished: updateAssocidatedGeneIds.");
        return namespaces;
    }

    private static void createOverlapout(HashMap<String, HashMap<String, GOclass>> namespaces, String root,
                                         File overlapout, int minsize, int maxsize) {
        System.out.println("Begin: createOverlapout (root: " + root + ")...");
        // conversion to array
        GOclass[] namespaceArr = namespaces.get(root).values().toArray(new GOclass[0]);
        ArrayList<String> result = new ArrayList<>();
        result.add("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");

        // int counter = 0;
        for (int i = 0; i < namespaceArr.length - 1; i++) {
            // jump GO's with minsize > size > maxsize (First Pair)
            if (namespaceArr[i].associatedGenes.size() >= minsize
                    && namespaceArr[i].associatedGenes.size() <= maxsize) {
                // counter++;
                // System.out.println("counter: " + counter);
                for (int j = i + 1; j < namespaceArr.length; j++) {
                    // jump GO's with minsize > size > maxsize (Second Pair)
                    if (namespaceArr[j].associatedGenes.size() >= minsize
                            && namespaceArr[j].associatedGenes.size() <= maxsize) {
                        // for each pair of GO's:
                        String term1 = namespaceArr[i].id;
                        String term2 = namespaceArr[j].id;
                        boolean is_relative = false;
                        double pathLengthBestResult = Integer.MAX_VALUE;
                        int num_overlapping = 0;
                        double max_ov_percent = 0;

                        // create intersection. Copy-constructor i & retain strings occurring in j
                        HashSet<String> intersection = new HashSet<String>(namespaceArr[i].associatedGenes);
                        intersection.retainAll(namespaceArr[j].associatedGenes);
                        num_overlapping = intersection.size();

                        if (num_overlapping > 0) {
                            double min_associatedGenes = (double) namespaceArr[i].associatedGenes.size();
                            if (namespaceArr[j].associatedGenes.size() < namespaceArr[i].associatedGenes.size()) {
                                min_associatedGenes = namespaceArr[j].associatedGenes.size();
                            }
                            max_ov_percent = (double) num_overlapping / min_associatedGenes;
                            max_ov_percent = (double) Math.round(max_ov_percent * 100 * 100d) / 100d;
                            // System.out.println("num_overlapping: " + num_overlapping + " i and j :"
                            // + namespaceArr[i].associatedGenes.size() + " and "
                            // + namespaceArr[j].associatedGenes.size());
                            // System.out.println("Rounded result: " + max_ov_percent);

                            // is_relative: true if an associated DAG entry to term1 is ascendent or
                            // descendant of the one associated to term2, false otherwise
                            if (namespaceArr[j].is_a_updated.contains(namespaceArr[i].id)) {
                                is_relative = true;
                                // calculate pathLength (is relative)
                                pathLengthBestResult = getPathBetweenRelatives(namespaces.get(root), namespaceArr[i].id,
                                        namespaceArr[j].id);
                            } else if (namespaceArr[i].is_a_updated.contains(namespaceArr[j].id)) {
                                is_relative = true;
                                // calculate pathLength (is relative)
                                pathLengthBestResult = getPathBetweenRelatives(namespaces.get(root), namespaceArr[j].id,
                                        namespaceArr[i].id);
                                pathLengthBestResult++;
                            } else {
                                // calculate pathLength (is NOT relative)
                                // Gather list containing all GO's occuring in i and j
                                ArrayList<String> commonGOs = new ArrayList<String>();
                                for (String is_a_item : namespaceArr[i].is_a_updated) {
                                    if (namespaceArr[j].is_a_updated.contains(is_a_item)) {
                                        commonGOs.add(is_a_item);
                                    }
                                }
                                if (!commonGOs.isEmpty()) {
                                    // Find the shortest path between all common GO's
                                    for (String is_a_item : commonGOs) {
                                        double currPathLength = getPathBetweenRelatives(namespaces.get(root), is_a_item,
                                                namespaceArr[i].id);
                                        currPathLength += getPathBetweenRelatives(namespaces.get(root), is_a_item,
                                                namespaceArr[j].id);
                                        currPathLength++;
                                        if (currPathLength < pathLengthBestResult) {
                                            pathLengthBestResult = currPathLength;
                                        }
                                    }
                                } else {
                                    // first common ancestor = root --> path_length = length to root
                                    pathLengthBestResult = (getLengthToRoot(namespaces.get(root), namespaceArr[j].id)
                                            + getLengthToRoot(namespaces.get(root), namespaceArr[i].id));
                                }
                            }

                            result.add(term1 + "\t" + term2 + "\t" + is_relative + "\t" + (int) pathLengthBestResult
                                    + "\t" + num_overlapping + "\t" + max_ov_percent);
                        }
                    }
                }
            }
        }
        // write result to its savespace
        System.out.println("Finished: createOverlapout");
        writeArrayList(result, overlapout);
    }

    public static double getPathBetweenRelatives(HashMap<String, GOclass> namespace, String parent, String child) {
        // System.out.println("****|getPathBetweenRelatives| GO_id: " + term1 + " and: " + term2 + "****");
        if (namespace.get(child).is_a_original.isEmpty()) {
            // This path leads not to a result
            return Integer.MAX_VALUE;
        }
        if (namespace.get(child).is_a_original.contains(parent)) {
            return 0;
        } else {
            double bestResult = Integer.MAX_VALUE;
            for (String go : namespace.get(child).is_a_original) {
                double tempResult = getPathBetweenRelatives(namespace, parent, go);
                if (tempResult < Integer.MAX_VALUE) {
                    if (tempResult < bestResult) {
                        bestResult = tempResult;
                    }
                }
            }
            return 1 + bestResult;
        }
    }

    public static HashMap<String, HashMap<String, GOclass>> parseAssociatedGeneIds(String mappingtype, File mapping, HashMap<String, HashMap<String, GOclass>> namespaces, String root) {
        // mappingtype defines how mapping is processed
        if (mappingtype.equals("go")) {
            try {
                System.out.println("Begin: parseAssociatedGeneIds (mappingtype: " + mappingtype + ")...");
                BufferedReader br = new BufferedReader(
                        new InputStreamReader(new GZIPInputStream(new FileInputStream(mapping)), "UTF-8"));
                String line = "";
                while ((line = br.readLine()) != null) {
                    if (!line.startsWith("!")) {
                        // ignore comments ("!")
                        String gene_id = line.split("\t")[2];
                        String modifier = line.split("\t")[3];
                        String GO_id = line.split("\t")[4];
                        // error checking (modifier is optional and can be empty)
                        if (gene_id.isEmpty() || GO_id.isEmpty()) {
                            System.err.println("|parseAssociatedGeneIds| Empty value found. Line: " + line);
                        } else if (modifier.isEmpty()) {
                            // use value only if modifier is empty (see assignment)
                            if (namespaces.get(root).containsKey(GO_id)) {
                                namespaces.get(root).get(GO_id).associatedGenes.add(gene_id);
                                // System.out.println("Gene_id: "+gene_id+" was added to GeneOntology "+GO_id);
                            } else {
                                // System.err.println(
                                // "|parseAssociatedGeneIds| GO_id " + GO_id + " was not found in root: " +
                                // root);
                            }
                        }
                        // System.out.println("gene_id: " + gene_id + " modifier: " + modifier + "
                        // GO_id: " + GO_id);
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else if (mappingtype.equals("ensembl")) {
            System.out.println("Begin: parseAssociatedGeneIds (mappingtype: " + mappingtype + ")...");
            try {
                BufferedReader br = new BufferedReader(new FileReader(mapping));
                br.readLine(); // this will read and skip the first line
                String line;
                while ((line = br.readLine()) != null) {
                    // System.out.println("Line: " + line);
                    String gene_id = line.split("\t")[1];
                    // System.out.println("Gene id: "+gene_id);
                    String[] GO_ids = line.split("\t")[2].split("\\|");
                    // System.out.println("GO_ids:" + Arrays.toString(GO_ids));
                    // error checking
                    if (gene_id.isEmpty()) {
//                        System.err.println("|parseAssociatedGeneIds| Empty gene_id found. Line: " + line);
                    } else if (GO_ids.length < 1) {
                        System.err.println("|parseAssociatedGeneIds| GO_ids is empty! Line: " + line);
                    } else {
                        for (String currGOid : GO_ids) {
                            if (namespaces.get(root).containsKey(currGOid)) {
                                namespaces.get(root).get(currGOid).associatedGenes.add(gene_id);
                                // System.out.println("Gene_id: " + gene_id + " was added to GeneOntology " +
                                // currGOid);
                            } else {
                                // System.err.println("|parseAssociatedGeneIds| currGOid " + currGOid
                                // + " was not found in root: " + root);
                            }
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            System.err
                    .println("|parseAssociatedGeneIds| Mappingtype must correspond to [go/ensembl]. Found mappingtype: "
                            + mappingtype);
        }
        System.out.println("Finished: parseAssociatedGeneIds.");
        return namespaces;
    }

    public static HashMap<String, HashMap<String, GOclass>> parseOboFile(File obo) {
        // Returns a hashmap with key = namespace, value = GOclasses belonging to the
        // respective namespace
        HashMap<String, HashMap<String, GOclass>> namespaces = new HashMap<String, HashMap<String, GOclass>>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(obo));
            String line = "";
            System.out.println("Begin: parseOboFile...");
            boolean termOpen = false;
            String id = "";
            String name = "";
            String namespace = "";
            HashSet<String> is_a = new HashSet<String>();
            boolean is_obsolete = false;
            while ((line = br.readLine()) != null) {
                // System.out.println("currLine: " + line);
                if (termOpen) {
                    // gather information of current term
                    // id: GO:0000001
                    // name: mitochondrion inheritance
                    // namespace: biological_process
                    // def: "The distribu...
                    // synonym: "mitochondrial inheritance" EXACT []
                    // is_a: GO:0048308 ! organelle inheritance
                    // is_a: GO:0048311 ! mitochondrion distribution
                    // evtl. is_obsolete: true
                    if (line.startsWith("id:")) {
                        id = line.split(" ")[1];
                        // error checking
                        if (line.split(" ").length > 2) {
                            System.err.println("More than one space found in " + line);
                        }
                    } else if (line.startsWith("name:")) {
                        name = line.substring(6);
                        // System.out.println("name: " + name);
                        // error checking
                        if (name.isEmpty()) {
                            System.err.println("Empty name found in line: " + line);
                        }
                    } else if (line.startsWith("namespace:")) {
                        namespace = line.split(" ")[1];
                        // error checking
                        if (line.split(" ").length > 2) {
                            System.err.println("More than one space found in " + line);
                        }
                    } else if (line.startsWith("is_a:")) {
                        is_a.add(line.split(" ")[1]);
                        // no error checking!
                    } else if (line.startsWith("is_obsolete: true")) {
                        is_obsolete = true;
                    } else if (line.equals("")) {
                        // End of current term --> Save term to namespaces
                        // System.out.println("Reached term end. namespace: " + namespace + " id: " + id
                        // + " is_a.size(): " + is_a.size());
                        if (!is_obsolete) {
                            if (!namespaces.containsKey(namespace)) {
                                namespaces.put(namespace, new HashMap<String, GOclass>());
                            }
                            // ADD all the "is_a-GOclasses" from current GOclass
                            // for (String parentId : is_a) {
                            // if (!namespaces.get(namespace).containsKey(parentId)) {
                            // namespaces.get(namespace).put(parentId, new GOclass(parentId, namespace));
                            // }
                            // namespaces.get(namespace).get(parentId).contains_GOs.add(id);
                            // }
                            // ADD the current GOclass itself to namespaces
                            if (!namespaces.get(namespace).containsKey(id)) {
                                namespaces.get(namespace).put(id, new GOclass(id, name, namespace));
                            }
                            // update List of is_a's
                            for (String parentId : is_a) {
                                namespaces.get(namespace).get(id).is_a_updated.add(parentId);
                                namespaces.get(namespace).get(id).is_a_original.add(parentId);
                            }
                        }
                        // Reset values to default
                        termOpen = false;
                        id = "";
                        namespace = "";
                        is_a = new HashSet<String>();
                        is_obsolete = false;
                    }
                } else {
                    // wait for new term to open
                    if (line.equals("[Term]")) {
                        termOpen = true;
                    }
                }
            }
            br.close();
        } catch (

                Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished: parseOboFile.");
        return namespaces;
    }

    public static void writeNamespaces(HashMap<String, HashMap<String, GOclass>> namespaces, String root, File
            output) {
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(output, false)));
            out.println("|printNamespaces|");
            // for (String ns : namespaces.keySet()) {
            for (String id_goclass : namespaces.get(root).keySet()) {
                out.print("root: " + root + " GOclass_id: " + id_goclass + "\tis_a:");
                for (String is_a_temp : namespaces.get(root).get(id_goclass).is_a_updated) {
                    out.print(" " + is_a_temp);
                }
                out.print(" AssGenes-size:" + namespaces.get(root).get(id_goclass).associatedGenes.size() + " (");
                for (String g : namespaces.get(root).get(id_goclass).associatedGenes) {
                    out.print(g + ", ");
                }
                out.print(")");
                out.println();
            }
            // }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void writeArrayList(ArrayList<String> result, File location) {
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(location, false)));
            for (String string : result) {
                out.println(string);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static HashMap<String, HashMap<String, GOclass>> updateGOlinks(
            HashMap<String, HashMap<String, GOclass>> namespaces, String root) {
        System.out.println("Begin: updateGOlinks (root: " + root + ")...");
        // returns an updated version of namespaces containing ALL links (not only is_a)
        for (String GO_id : namespaces.get(root).keySet()) {
            // System.out.println("Updating GO_id: " + GO_id);
            // for any GOclass in namespaces, call getParentalIsALinks
            namespaces.get(root).get(GO_id).is_a_updated.addAll(getParentalIsALinks(namespaces.get(root), GO_id));
        }
        System.out.println("Finished: updateGOlinks.");
        return namespaces;
    }

    public static int getLengthToRoot(HashMap<String, GOclass> GOclassesInNamespace, String GO_id) {
        // System.out.println("****|getLengthToRoot| GO_id: " + GO_id + "****");
        if (GOclassesInNamespace.get(GO_id).level > -1) {
            return GOclassesInNamespace.get(GO_id).level;
        } else {
            if (GOclassesInNamespace.get(GO_id).is_a_original.size() == 0) {
                // System.out.println("is_a_original SIZE = 0");
                GOclassesInNamespace.get(GO_id).level = 1;
                return 1;
            } else {
                int bestResult = Integer.MAX_VALUE;
                for (String go : GOclassesInNamespace.get(GO_id).is_a_original) {
                    // System.out.println(GO_id + " contains " + go);
                    // System.out.println("is_a_original.size: " +
                    // GOclassesInNamespace.get(GO_id).is_a_original.size());
                    int tempResult = getLengthToRoot(GOclassesInNamespace, go);
                    GOclassesInNamespace.get(go).level = tempResult;
                    if (tempResult < bestResult) {
                        bestResult = tempResult;
                    }
                }
                return 1 + bestResult;
            }
        }
    }

    public static HashSet<String> getParentalIsALinks(HashMap<String, GOclass> GOclassesInNamespace, String
            GO_id) {
        // RECURSIVE
        // result = all parental is_a links of respective GO_id
        HashSet<String> result = new HashSet<String>();
        if (GOclassesInNamespace.get(GO_id).is_a_updated.size() != 0) {
            // parents exist
            for (String parentId : GOclassesInNamespace.get(GO_id).is_a_updated) {
                result.addAll(GOclassesInNamespace.get(parentId).is_a_updated);
                result.addAll(getParentalIsALinks(GOclassesInNamespace, parentId));
            }
        } else {
            // parents don't exist
        }
        return result;
    }

}
