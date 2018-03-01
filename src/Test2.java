import java.util.ArrayList;

public class Test2 {


    public static void main(String[] args) {
        System.out.println("start");

        ArrayList<String> test = new ArrayList<>();

        test.add("eins");
        test.add("zwei");

        System.out.println(test.toString());

        try {
            System.out.println("0: "+test(0));
        } catch (Exception e) {
//            e.printStackTrace();
        }

        try {
            System.out.println("1: "+test(1));
        } catch (Exception e) {
//            e.printStackTrace();
        }

    }


    public static int test(int zahl) throws Exception {
        // error wenn zahl != 1
        int result = 0;
        if (zahl == 1) {
            result = 23;
        } else {
            throw new Exception("This path leads not to a result.");
        }
        return result;
    }

}
