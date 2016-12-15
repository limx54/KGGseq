package org.cobi.bayes;

import cern.colt.list.DoubleArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import org.apache.log4j.Logger;
<<<<<<< HEAD
import org.cobi.kggseq.Constants;
=======
>>>>>>> origin/master

import org.cobi.util.file.LocalFileFunc;

public class Bayes implements Constants {

    String resourcePath;
    private static final Logger LOG = Logger.getLogger(Bayes.class);

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub
/*
        The population mean score for imputation of missing value.
        GWAVA_Region	0.153309
GWAVA_TSS	0.205698
GWAVA_Unmatched	0.33259
CADD_CScore	-0.0784622
DANN	0.291601
FATHMM-MKL	0.00403947
FunSeq	1
FunSeq2	0.287838
GWAS3D	2.80186
SuRFR	9.44965        
         */

<<<<<<< HEAD
        double[] annotationScoreDouble = new double[]{0.153309, 0.205698, 0.33259, -0.0784622, 0.291601, 0.00403947, 1, -0.540851872, 2.80186, 9.44965};

        Bayes bayesScoreGenerator = new Bayes("resources", "hg19");

        bayesScoreGenerator.readResource(false);
        double causalscore, neutralscore;
        for (int i = 0; i < annotationScoreDouble.length; i++) {
            causalscore = bayesScoreGenerator.getCausalscore(annotationScoreDouble[i], i);
            neutralscore = bayesScoreGenerator.getNeutralscore(annotationScoreDouble[i], i);
            System.out.println(causalscore + "\t" + neutralscore + "\t" + (causalscore / neutralscore) + "\t" + (causalscore / (causalscore + neutralscore)));
=======
        double[] annotationScoreDouble = new double[]{-0.111646, 0.1806, 0.381241, 0.81, 0.6, 0.37, 17.519};
        double[] annotationScore = new double[annotationScoreDouble.length];
        for (int i = 0; i < annotationScore.length; i++) {
            annotationScore[i] = (Double) annotationScoreDouble[i];
>>>>>>> origin/master
        }
        // double cell_p = bayesScoreGenerator.getCellSpecificScore("6", 127663116);
        //System.out.println(bayesScoreGenerator.getBayesScore(annotationScore, cell_p));
    }

    public void readResource(boolean cellTypeSpec) {
        causalScores = getDistributionProbability(resourcePath + "/" + causalFilePath + "/");
        neutralScores = getDistributionProbability(resourcePath + "/" + neutralFilePath + "/");
        if (cellTypeSpec) {
            cellTypeScores = getCellSpecificityElements();
        }
    }

    public String causalFilePath = "all_causal_distribution";
    public String neutralFilePath = "all_neutral_distribution";
    public String cellLineName = "GM12878";
    public String cellEncodeInfoName[] = {"H3K4me1", "H3K36me3", "DNase", "H3K79me2", "H3K9me3", "H3K27me3", "H3K4me2", "H3K4me3", "H3K36me3"};

    public ArrayList<double[]> causalScores = new ArrayList<double[]>();
    public ArrayList<double[]> neutralScores = new ArrayList<double[]>();
    public ArrayList<HashMap<String, ArrayList<double[]>>> cellTypeScores;
<<<<<<< HEAD

=======
>>>>>>> origin/master
    public StringBuffer sb;
    public String refGenomeVersion;

    public Bayes(String resourcePath, String refGenomeVersion) {
        this.resourcePath = resourcePath;
        this.refGenomeVersion = refGenomeVersion;
    }

    /*
    public void changeFeatureNum(String[] inputList) {
        if (inputList != null) {
            featureNum = new String[inputList.length];
            for (int i = 0; i < inputList.length; i++) {
                featureNum[i] = inputList[i];
            }
        }
    }
     */
    public void changeCellLineName(String inputName) {
        if (inputName != null) {
            cellLineName = inputName;
        }
    }

    public float getCellSpecificScore(String ChromeIndex, int posIndex) {
        float score = 0;
        double List_Hit[] = new double[6];
        double H3K79me2_centrality = -1;
        for (int ii = 0; ii < List_Hit.length; ii++) {
            if (cellTypeScores.get(ii).containsKey(ChromeIndex)) {
                double[] annotationPosition = new double[2 * cellTypeScores.get(ii).get(ChromeIndex).size()];
                for (int i = 0; i < cellTypeScores.get(ii).get(ChromeIndex).size(); i++) {
                    annotationPosition[2 * i] = cellTypeScores.get(ii).get(ChromeIndex).get(i)[0];
                    annotationPosition[2 * i + 1] = cellTypeScores.get(ii).get(ChromeIndex).get(i)[1];
                }
                int searchIndex = getBinarySearchRegion(annotationPosition, (float) posIndex, 0, annotationPosition.length / 2);
                if (searchIndex > 0) {
                    List_Hit[ii] = 1;
                    if (ii == 3) {
                        H3K79me2_centrality = Math.abs(cellTypeScores.get(ii).get(ChromeIndex).get(searchIndex)[3]
                                - (posIndex - cellTypeScores.get(ii).get(ChromeIndex).get(searchIndex)[0]));
                    }
                } else {
                    List_Hit[ii] = 0;
                }
            } else {
                List_Hit[ii] = 0;
            }
        }
        double H3K4me1_hit = List_Hit[0];
        double H3K36me3_hit = List_Hit[1];
        double DNase_hit = List_Hit[2];
        double H3K79me2_hit = List_Hit[3];
        double H3K9me3_hit = List_Hit[4];
        double H3K27me3_hit = List_Hit[5];

        double[] List_Score = new double[3];
        for (int a = 0; a < List_Score.length; a++) {
            int ii = a + 6;
            if (cellTypeScores.get(ii).containsKey(ChromeIndex)) {
                double[] annotationPosition = new double[2 * cellTypeScores.get(ii).get(ChromeIndex).size()];
                for (int i = 0; i < cellTypeScores.get(ii).get(ChromeIndex).size(); i++) {
                    annotationPosition[2 * i] = cellTypeScores.get(ii).get(ChromeIndex).get(i)[0];
                    annotationPosition[2 * i + 1] = cellTypeScores.get(ii).get(ChromeIndex).get(i)[1];
                }
                int searchIndex = getBinarySearchRegion(annotationPosition, (float) posIndex, 0, annotationPosition.length / 2);
                if (searchIndex > 0) {
                    List_Score[a] = cellTypeScores.get(ii).get(ChromeIndex).get(searchIndex)[2];
                } else {
                    List_Score[a] = 0;
                }
            } else {
                List_Score[a] = 0;
            }
        }
        double H3K4me2_score = List_Score[0];
        double H3K4me3_score = List_Score[1];
        double H3K36me3_score = List_Score[2];
        score = (float) (1 / (1 + Math
                .exp(-(-0.5339527052 + 1.0513562209 * H3K4me1_hit + 1.5659681399 * H3K36me3_hit + 1.2131942069 * DNase_hit + 0.9750312605 * H3K79me2_hit + -0.4843821400 * H3K9me3_hit
                        + 1.5150317212 * H3K27me3_hit + 0.0008691201 * H3K4me2_score + 0.0003089830 * H3K4me3_score + 0.0043517819 * H3K36me3_score + -0.0001497833 * H3K79me2_centrality))));
        if (score < 0.3696304) {
            score = (float) 0.3696304;
        }
        return score;
    }

    public int getBinarySearchRegion(double[] List, float posIndex, int start, int end) {
        int mid = start + (end - start) / 2;
        if (end < start || posIndex > List[List.length - 1] || posIndex < List[0]) {
            return -1;
        }
        if (List[2 * mid] > posIndex) {
            return getBinarySearchRegion(List, posIndex, start, mid - 1);
        } else if (List[2 * mid + 1] < posIndex) {
            return getBinarySearchRegion(List, posIndex, mid + 1, end);
        } else {
            return mid;
        }
    }

<<<<<<< HEAD
    public double getCausalscore(double s, int index) {
        double causalscore = causalScores.get(2 * index + 1)[getBinarySearch(causalScores.get(2 * index), s, 0, causalScores.get(2 * index).length - 1)];
        return causalscore;
    }

    public double getNeutralscore(double s, int index) {
        double causalscore = neutralScores.get(2 * index + 1)[getBinarySearch(neutralScores.get(2 * index), s, 0, neutralScores.get(2 * index).length - 1)];
        return causalscore;
    }

    public double[] getBayesScoreCompsit(float[] annotationScore, int[] effecIndex, double[] baye1NonCodingPredic, double[] baye2NonCodingPredic) {
=======
    public String getBayesScore(double[] annotationScore, double cell_p) {
>>>>>>> origin/master
        double composite_p = 1;
        double bfFactor = 1;
        double[] result = new double[2];
        Arrays.fill(result, Double.NaN);
        double causalscore;
        double neutralscore;
        int index;

        //StringBuilder sbTest1 = new StringBuilder();
        //StringBuilder sbTest2 = new StringBuilder();
        for (int i = 0; i < effecIndex.length; i++) {
            causalscore = 1;
            neutralscore = 1;
            index = effecIndex[i];
            if (Double.isNaN(annotationScore[index])) {
                //tmpScore = baye1NonCodingPredic[effecIndex[i]];
                //causalscore = causalScores.get(2 * i + 1)[getBinarySearch(causalScores.get(2 * i), tmpScore, 0, causalScores.get(2 * i).length - 1)];
                //neutralscore = neutralScores.get(2 * i + 1)[getBinarySearch(neutralScores.get(2 * i), tmpScore, 0, neutralScores.get(2 * i).length - 1)];
                bfFactor *= baye1NonCodingPredic[index];
                composite_p *= baye2NonCodingPredic[index];
            } else {
<<<<<<< HEAD
                causalscore = causalScores.get(2 * index + 1)[getBinarySearch(causalScores.get(2 * index), annotationScore[index], 0, causalScores.get(2 * index).length - 1)];
                neutralscore = neutralScores.get(2 * index + 1)[getBinarySearch(neutralScores.get(2 * index), annotationScore[index], 0, neutralScores.get(2 * index).length - 1)];
                bfFactor *= (causalscore / neutralscore);
                composite_p *= (causalscore / (causalscore + neutralscore));
=======
                causalscore = causalScores.get(2 * i + 1)[getBinarySearch(causalScores.get(2 * i), annotationScore[i], 0, causalScores.get(2 * i).length - 1)];
                neutralscore = neutralScores.get(2 * i + 1)[getBinarySearch(neutralScores.get(2 * i), annotationScore[i], 0, neutralScores.get(2 * i).length - 1)];
//				sb.append(annotationScore[t] + "\t" + causalscore + "|" + neutralscore + "\t");
                sb.append(annotationScore[i] + "\t");
>>>>>>> origin/master
            }
            // System.out.println(i+" "+(causalscore / neutralscore)+" "+(causalscore / (causalscore + neutralscore)));
        }

        result[0] = bfFactor;
        result[1] = composite_p;
        return result;

    }

    public ArrayList<HashMap<String, ArrayList<double[]>>> getCellSpecificityElements() {
        ArrayList<HashMap<String, ArrayList<double[]>>> cellSpecificScore = new ArrayList<HashMap<String, ArrayList<double[]>>>();
        for (int i = 0; i < cellEncodeInfoName.length; i++) {
            File rsFile = new File(resourcePath + refGenomeVersion + "/all_cell_signal/" + cellLineName + "-" + cellEncodeInfoName[i] + ".narrowPeak.sorted.gz");
            //File rsFile = new File("C:\\Users\\mulin0424\\Desktop\\PRVCS\\resources\\" + genomeVersion + "\\all_cell_signal\\" + cellLineName + "-" + cellEncodeInfoName[i] + ".narrowPeak.sorted.gz");
            try {
                HashMap<String, ArrayList<double[]>> thisHashMap = new HashMap<String, ArrayList<double[]>>();
                BufferedReader br = LocalFileFunc.getBufferedReader(rsFile.getCanonicalPath());
                String currentLine;
                while ((currentLine = br.readLine()) != null) {
                    if (currentLine.trim().length() == 0) {
                        continue;
                    }
                    String sp[] = currentLine.trim().split("\t");
                    String key = "";
                    if (sp[0].startsWith("chr") || sp[0].startsWith("Chr")) {
                        key = sp[0].trim().substring(3, sp[0].trim().length());
                    } else {
                        key = sp[0].trim();
                    }
                    double[] element = new double[]{Double.parseDouble(sp[1].trim()), Double.parseDouble(sp[2].trim()), Double.parseDouble(sp[4].trim()),
                        Double.parseDouble(sp[9].trim())};
                    if (thisHashMap.containsKey(key)) {
                        thisHashMap.get(key).add(element);
                    } else {
                        ArrayList<double[]> thisScore = new ArrayList<double[]>();
                        thisScore.add(element);
                        thisHashMap.put(key, thisScore);
                    }

                }
                cellSpecificScore.add(thisHashMap);
                br.close();
            } catch (Exception e) {
                // TODO Auto-generated catch block
                LOG.error(rsFile.toString() + " doesn't exist!");
            }

        }
        return cellSpecificScore;
    }

    public int getBinarySearch(final double[] List, double key, int start, int end) {
        if (key >= List[List.length - 1]) {
            return List.length - 1;
        } else if (key <= List[0]) {
            return 0;
        }
        if (end < start) {
            double mid = (List[end] + List[start]) / 2;
            if (key >= mid) {
                return start;
            }
            return end;
        }
        int mid = start + (end - start) / 2;
        if (key > List[mid]) {
            return getBinarySearch(List, key, mid + 1, end);
        } else if (key < List[mid]) {
            return getBinarySearch(List, key, start, mid - 1);
        }
        return mid;

    }

<<<<<<< HEAD
    public ArrayList<double[]> getDistributionProbability(String filePath) {
        ArrayList<double[]> scores = new ArrayList<double[]>();
        DoubleArrayList thisScore = new DoubleArrayList();
        DoubleArrayList thisProbability = new DoubleArrayList();
        File resourceFile = null;
        for (int i = 0; i < ReguPredicNames.length; i++) {
=======
    public ArrayList<double[]> getDistributionProbability(String FilePath) {
        ArrayList<double[]> scores = new ArrayList<double[]>();
        DoubleArrayList thisScore = new DoubleArrayList();
        DoubleArrayList thisProbability = new DoubleArrayList();
        for (int i = 0; i < featureNum.length; i++) {
>>>>>>> origin/master
            thisScore.clear();
            thisProbability.clear();
            try {
                resourceFile = new File(filePath + ReguPredicNames[i] + ".dis");
                if (!resourceFile.exists()) {
                    LOG.error(resourceFile.getCanonicalPath() + " does not exist!");
                    continue;
                }
                BufferedReader br = new BufferedReader(new FileReader(resourceFile));
                while (br.ready()) {
                    String str = br.readLine();
                    String[] splitStr = str.trim().split(" ");
                    if (splitStr.length == 3) {
                        /*
                        if (ReguPredicNames[i].equals("FunSeq2_score")) {
                            thisScore.add(Math.pow(10, Double.parseDouble(splitStr[0])));
                        } else {
                            thisScore.add((double) Double.parseDouble(splitStr[0]));
                        }
                         */
                        thisScore.add((double) Double.parseDouble(splitStr[0]));
                        thisProbability.add((double) Double.parseDouble(splitStr[2]));
                    }
                }
                double[] floatThisScore = new double[thisScore.size()];
                double[] floatThisProability = new double[thisProbability.size()];
                for (int t = 0; t < floatThisScore.length; t++) {
                    floatThisScore[t] = thisScore.getQuick(t);
                }
                for (int t = 0; t < floatThisScore.length; t++) {
                    floatThisProability[t] = thisProbability.getQuick(t);
                }

                scores.add(floatThisScore);
                scores.add(floatThisProability);
                br.close();
            } catch (Exception e) {
                LOG.error(e);
            }
        }
        return scores;
    }

}
