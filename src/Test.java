import java.io.*;
import java.util.HashMap;

public class Test {

    public static void main(String[] args) {

        File mysolution = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt5\\GOEnrichment\\result\\out.tsv");
        File theirsolution = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt5\\GOEnrichment\\simul_exp_go_bp_ensembl_min50_max500.enrich.out");

        HashMap<String, String> myParse = parseFile(mysolution);
        HashMap<String, String> theirParse = parseFile(theirsolution);
        System.out.println("myParse.size: " + myParse.size());
        System.out.println("theirParse.size: " + theirParse.size());

        int counter = 0;
        for (String myGeneId : myParse.keySet()) {
            if (theirParse.containsKey(myGeneId)) {
                // my geneId is contained in their solution
                int myPipe = myParse.get(myGeneId).length() - myParse.get(myGeneId).replace("|", "").length();
                int theirPipe = theirParse.get(myGeneId).length() - theirParse.get(myGeneId).replace("|", "").length();

                if (myPipe != theirPipe) {
                    counter++;
                    System.err.println("Amount of pipes not equal. (Nr.:" + counter + ")");
                    System.err.println("Myyyy path (" + myPipe + "): " + myGeneId + " " + myParse.get(myGeneId));
                    System.err.println("Their path (" + theirPipe + "): " + myGeneId + " " + theirParse.get(myGeneId));
                }

            } else {
                System.err.println("Key " + myGeneId + " not found in their solution.");
            }
        }


    }

    public static HashMap<String, String> parseFile(File file) {
        System.out.println("Begin: parseFile..." + file.getAbsolutePath());
        HashMap<String, String> result = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            while ((line = br.readLine()) != null) {
                String[] lineArr = line.split("\\t");
                if (lineArr.length != 12) {
                    String geneId = lineArr[0];
                    String path = lineArr[12];
                    result.put(geneId, path);
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }
}
