import java.io.*;
import java.util.Arrays;
import java.util.List;

public class FragmentsExportMainClass {

    public static void main(String[] args) {
        File resultsFolder = new File(args[0]);
        int threadsNumber  = Integer.parseInt(args[1]);
        String ionsTypes   = args[2];
        try {
            new FragmentsExportMainClass(resultsFolder, threadsNumber, ionsTypes);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public FragmentsExportMainClass(File resultsFolder, int threadsNumber, String ionsTypes) throws IOException {
        List<String> ionsTypeArray = Arrays.asList(ionsTypes.split("_"));
        new ExportFragments(resultsFolder, threadsNumber, ionsTypeArray);
    }
}
