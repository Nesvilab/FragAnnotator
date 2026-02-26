import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.identification.identification_parameters.PtmSettings;
import com.compomics.util.experiment.identification.identification_parameters.SearchParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationSettings;

import java.io.*;
import java.util.ArrayList;

public class FragmentsExportMainClass {

    /**
     * Annotation setting
     */
    public AnnotationSettings annotationSettings = new AnnotationSettings();

    public static void main(String[] args) {

        File resultsFolder = new File(args[0]);
        int threadsNumber = Integer.parseInt(args[1]);

        try {
            new FragmentsExportMainClass(resultsFolder, threadsNumber);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public FragmentsExportMainClass(File resultsFolder, int threadsNumber) throws IOException {

        PTMFactory ptmFactory = PTMFactory.getInstance();
        SearchParameters searchParameters = new SearchParameters();

        ArrayList<String> modification = ptmFactory.getPTMs();
        PtmSettings ptmSettings = new PtmSettings();

        for(String fixedModification:modification){
            ptmSettings.addFixedModification(ptmFactory.getPTM(fixedModification));
        }

        for(String variableModification:modification){
            ptmSettings.addVariableModification(ptmFactory.getPTM(variableModification));
        }

        searchParameters.setPtmSettings(ptmSettings);

        searchParameters.setFragmentAccuracyType(SearchParameters.MassAccuracyType.PPM);

        searchParameters.setFragmentIonAccuracy(20.0);

        annotationSettings.setPreferencesFromSearchParameters(searchParameters);

        ExportFragments exportFragments = new ExportFragments(resultsFolder, annotationSettings, threadsNumber);

    }

}
