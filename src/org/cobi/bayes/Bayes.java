package org.cobi.bayes;

import cern.colt.list.DoubleArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import org.apache.log4j.Logger;

import org.cobi.util.file.LocalFileFunc;

public class Bayes {

    String genomeVersion;
    private static final Logger LOG = Logger.getLogger(Bayes.class);

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub

        double[] annotationScoreDouble = new double[]{-0.111646, 0.1806, 0.381241, 0.81, 0.6, 0.37, 17.519};
        double[] annotationScore = new double[annotationScoreDouble.length];
        for (int i = 0; i < annotationScore.length; i++) {
            annotationScore[i] = (Double) annotationScoreDouble[i];
        }
        Bayes bayesScoreGenerator = new Bayes("hg19");
        bayesScoreGenerator.featureNum = new String[]{"CADD_cscore", "DANN_score", "FunSeq_score", "FunSeq2_score", "GWAS3D_score", "GWAVA_TSS_score", "SuRFR_score",
            "Fathmm_MKL_score"};
        bayesScoreGenerator.readResource();
        double cell_p = bayesScoreGenerator.getCellSpecificScore("6", 127663116);
        System.out.println(bayesScoreGenerator.getBayesScore(annotationScore, cell_p));
    }

    public void readResource() {
        causalScores = getDistributionProbability(causalFilePath);
        neutralScores = getDistributionProbability(neutralFilePath);
        cellTypeScores = getCellSpecificityElements();
    }

    public String causalFilePath = "./resources/all_causal_distribution/";
    public String neutralFilePath = "./resources/all_neutral_distribution/";
    public String featureNum[] = {"CADD_cscore", "DANN_score", "FunSeq_score", "FunSeq2_score", "GWAS3D_score", "GWAVA_TSS_score", "SuRFR_score", "Fathmm_MKL_score"};
    public String cellLineName = "GM12878";
    public String cellEncodeInfoName[] = {"H3K4me1", "H3K36me3", "DNase", "H3K79me2", "H3K9me3", "H3K27me3", "H3K4me2", "H3K4me3", "H3K36me3"};

    public ArrayList<double[]> causalScores = new ArrayList<double[]>();
    public ArrayList<double[]> neutralScores = new ArrayList<double[]>();
    public ArrayList<HashMap<String, ArrayList<double[]>>> cellTypeScores;
    public StringBuffer sb;

    public Bayes(String genomeVersion) {
        this.genomeVersion = genomeVersion;
    }

    public void changeFeatureNum(String[] inputList) {
        if (inputList != null) {
            featureNum = new String[inputList.length];
            for (int i = 0; i < inputList.length; i++) {
                featureNum[i] = inputList[i];
            }
        }
    }

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
        Double H3K4me1_hit = List_Hit[0];
        Double H3K36me3_hit = List_Hit[1];
        Double DNase_hit = List_Hit[2];
        Double H3K79me2_hit = List_Hit[3];
        Double H3K9me3_hit = List_Hit[4];
        Double H3K27me3_hit = List_Hit[5];

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
        Double H3K4me2_score = List_Score[0];
        Double H3K4me3_score = List_Score[1];
        Double H3K36me3_score = List_Score[2];
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

    public String getBayesScore(double[] annotationScore, double cell_p) {
        double composite_p = 1;
        double bfFactor = 1;
        sb = new StringBuffer();
        for (int i = 0; i < annotationScore.length; i++) {
            double causalscore = 1;
            double neutralscore = 1;
            if (Double.isNaN(annotationScore[i])) {
                causalscore = Math.log(0.1);
                neutralscore = Math.log(0.9);
                sb.append(".\t");
//				sb.append("NaN" + "\t" + 0.1 + "|" + 0.9 + "\t");
            } else {
                causalscore = causalScores.get(2 * i + 1)[getBinarySearch(causalScores.get(2 * i), annotationScore[i], 0, causalScores.get(2 * i).length - 1)];
                neutralscore = neutralScores.get(2 * i + 1)[getBinarySearch(neutralScores.get(2 * i), annotationScore[i], 0, neutralScores.get(2 * i).length - 1)];
//				sb.append(annotationScore[t] + "\t" + causalscore + "|" + neutralscore + "\t");
                sb.append(annotationScore[i] + "\t");
            }
            bfFactor *= causalscore / neutralscore;
            composite_p *= causalscore / (causalscore + neutralscore);
        }
        if (Double.isNaN(cell_p)) {
            cell_p = 0.3696304;
        }
        sb.append(bfFactor + "\t" + composite_p + "\t");
        sb.append(cell_p + "\t" + composite_p * cell_p / 0.5);
        return sb.toString();

    }

    public ArrayList<HashMap<String, ArrayList<double[]>>> getCellSpecificityElements() {
        ArrayList<HashMap<String, ArrayList<double[]>>> cellSpecificScore = new ArrayList<HashMap<String, ArrayList<double[]>>>();
        for (int i = 0; i < cellEncodeInfoName.length; i++) {
            File rsFile = new File("./resources/" + genomeVersion + "/all_cell_signal/" + cellLineName + "-" + cellEncodeInfoName[i] + ".narrowPeak.sorted.gz");
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

    public ArrayList<double[]> getDistributionProbability(String FilePath) {
        ArrayList<double[]> scores = new ArrayList<double[]>();
        DoubleArrayList thisScore = new DoubleArrayList();
        DoubleArrayList thisProbability = new DoubleArrayList();
        for (int i = 0; i < featureNum.length; i++) {
            thisScore.clear();
            thisProbability.clear();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(FilePath + featureNum[i] + ".dis")));
                while (br.ready()) {
                    String str = br.readLine();
                    String[] splitStr = str.trim().split(" ");
                    if (splitStr.length == 3) {
                        if (featureNum[i].equals("FunSeq2_score")) {
                            thisScore.add(Math.pow(10, Double.parseDouble(splitStr[0])));
                        } else {
                            thisScore.add((double) Double.parseDouble(splitStr[0]));
                        }
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
                // TODO: handle exception
                System.out.println("File does not exist");
            }
        }
        return scores;
    }

}
