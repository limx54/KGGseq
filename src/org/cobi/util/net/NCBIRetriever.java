/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.net;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.SocketException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import javax.net.ssl.HttpsURLConnection;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.apache.http.Header;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import org.cobi.util.download.stable.HttpClient4DownloadTask;
import org.cobi.util.text.Util;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;

/**
 *
 * @author mxli
 */
public class NCBIRetriever {

    DocumentBuilderFactory domfac = DocumentBuilderFactory.newInstance();

    //example
    // http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[mesh]%29+AND+TGM6[tiab]&datetype=edat&retmax=100
    //http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[All%20Fields]%29+AND+TGM6[All%20Fields]&datetype=edat&retmax=100
    public String pubMedIDESearch(List<String> diseaseNames, List<String> keyValues, List<String> filterValues) throws Exception {
        String line = null;
        if (diseaseNames.isEmpty() || keyValues.isEmpty()) {
            return null;
        }
        StringBuilder content1 = new StringBuilder();
        content1.append(URLEncoder.encode("db", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("pubmed", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("term", "UTF-8"));
        content1.append("=");
        StringBuilder idList = new StringBuilder();
        int diseaseNameNum = diseaseNames.size();
        int keyValuesNum = keyValues.size();
        StringBuilder diseaseStr = new StringBuilder();
        for (int i = 0; i < diseaseNameNum; i++) {
            String diseaseName = diseaseNames.get(i);
            diseaseName = diseaseNames.get(i);
            if (diseaseName == null || diseaseName.trim().length() == 0) {
                continue;
            }
            if (diseaseName.indexOf(' ') >= 0) {
                //()
                diseaseStr.append("(");
                diseaseStr.append(diseaseName.replace(' ', '+'));
                diseaseStr.append("[tiab])");
            } else {
                diseaseStr.append(diseaseName);
                diseaseStr.append("[tiab]");
            }
            diseaseStr.append("+OR+");
        }

        StringBuilder geneStr = new StringBuilder();
        for (int i = 0; i < keyValuesNum; i++) {
            String geneName = keyValues.get(i);
            if (geneName == null || geneName.trim().length() == 0) {
                continue;
            }
            if (geneName.indexOf(' ') >= 0) {
                //()
                geneStr.append("(");
                geneStr.append(geneName.replace(' ', '+'));
                geneStr.append("[tiab])");
            } else {
                geneStr.append(geneName);
                geneStr.append("[tiab]");
            }
            geneStr.append("+OR+");
        }
        StringBuilder filterStr = new StringBuilder();
        if (filterValues != null) {
            int filterValuesNum = filterValues.size();
            for (int i = 0; i < filterValuesNum; i++) {
                String filterName = filterValues.get(i);
                if (filterName == null || filterName.trim().length() == 0) {
                    continue;
                }
                if (filterName.indexOf(' ') >= 0) {
                    //()
                    filterStr.append("(");
                    filterStr.append(filterName.replace(' ', '+'));
                    filterStr.append("[tiab])");
                } else {
                    filterStr.append(filterName);
                    filterStr.append("[tiab]");
                }
                filterStr.append("+OR+");
            }
        }


        content1.append("((");
        content1.append(diseaseStr.substring(0, diseaseStr.length() - 4));
        content1.append(")+AND+(");
        content1.append(geneStr.substring(0, geneStr.length() - 4));
        if (filterStr.length() > 0) {
            content1.append(")+AND+(");
            content1.append(filterStr.substring(0, filterStr.length() - 4));
        }
        content1.append("))");


        content1.append("&");
        content1.append(URLEncoder.encode("datetype", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("edat", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("retmax", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("100", "UTF-8"));

        // System.out.println(content1);
        try {
              TimeUnit.SECONDS.sleep(4);
              
             
            // URL of CGI-Bin script.
            URL url = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?");

            //http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[mesh]%29+AND+TGM6[tiab]&datetype=edat&retmax=100
          

            // URL connection channel.
            HttpsURLConnection urlConn = (HttpsURLConnection)url.openConnection();
            // Let the run-time system (RTS) know that we want input.
            urlConn.setDoInput(true);
            // Let the RTS know that we want to do output.
            urlConn.setDoOutput(true);
            // No caching, we want the real thing.
            //urlConn.setUseCaches(false);
            //set timeout time
            urlConn.setConnectTimeout(700000);

            //System.out.println(content1);
            try ( //send data
            // Send quest output.
                    OutputStreamWriter wr = new OutputStreamWriter(urlConn.getOutputStream())) {
                //System.out.println(content1);
                wr.write(content1.toString());
                wr.flush();
            }
 
               
            /*
            BufferedReader rd = new BufferedReader(new InputStreamReader(urlConn.getInputStream()));
            while ((line = rd.readLine()) != null) {
            line = line.trim();
            if (line.length() == 0) {
            continue;
            }
            System.out.println(line);
            }
            rd.close();
             * 
             */


            DocumentBuilder dombuilder = domfac.newDocumentBuilder();
            Document doc = dombuilder.parse(urlConn.getInputStream());
            Element root = doc.getDocumentElement();

            NodeList nodes = root.getElementsByTagName("Count");
            if (nodes.getLength() == 1) {
                Element e = (Element) nodes.item(0);
                Text t = (Text) e.getFirstChild();
                if (t == null) {
                    return "";
                } else {
                    int count = Util.parseInt(t.getNodeValue());
                    if (count == 0) {
                        return "";
                    }
                }
            }
            /*
            nodes = root.getElementsByTagName("IdList");
            if (nodes.getLength() == 1) {
            Node item = nodes.item(0);
            for (Node node = item.getFirstChild(); node != null; node = node.getNextSibling()) {
            if (node.getNodeType() == Node.ELEMENT_NODE) {
            if (node.getNodeName().equals("Id")) {
            String name = node.getFirstChild().getNodeValue();
            System.out.println(name);
            }
            }
            }
            }
             * 
             */


            nodes = root.getElementsByTagName("Id");
            for (int i = 0; i < nodes.getLength(); i++) {
                Element subNode = (Element) nodes.item(i);
                Text t = (Text) subNode.getFirstChild();
                idList.append(t.getNodeValue());
                idList.append(';');
            }

            // System.out.println(idList.toString());
            if (idList.length() > 0) {
                return idList.substring(0, idList.length() - 1);
            } else {
                return "";
            }

        } catch (SocketException ex) {
            System.out.print("... failed!");
            return null;
        } catch (IOException ex) {
            System.out.print("... failed!");
            //ex.printStackTrace();
            return null;
        }

    }
    
         //example
    // http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[mesh]%29+AND+TGM6[tiab]&datetype=edat&retmax=100
    //http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[All%20Fields]%29+AND+TGM6[All%20Fields]&datetype=edat&retmax=100
    public String pubMedIDESearch(List<String> diseaseNames, String keyValue) throws Exception {
        String line = null;

        int index1 = 0;
        int index2 = 0;

        StringBuilder content1 = new StringBuilder();
        content1.append(URLEncoder.encode("db", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("pubmed", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("term", "UTF-8"));
        content1.append("=");
        StringBuilder idList = new StringBuilder();
        int diseaseNameNum = diseaseNames.size();

        if (diseaseNameNum > 1) {
            String diseaseName = diseaseNames.get(0);
            content1.append("%28");
            if (diseaseName.indexOf(' ') >= 0) {
                //()
                content1.append("%28");
                content1.append(diseaseName.replace(' ', '+'));
                content1.append("[tiab]%29");
            } else {
                content1.append(diseaseName);
                content1.append("[tiab]");
            }
            for (int i = 1; i < diseaseNameNum; i++) {
                content1.append("+OR+");
                diseaseName = diseaseNames.get(i);
                if (diseaseName.indexOf(' ') >= 0) {
                    //()
                    content1.append("%28");
                    content1.append(diseaseName.replace(' ', '+'));
                    content1.append("[tiab]%29");
                } else {
                    content1.append(diseaseName);
                    content1.append("[tiab]");
                }
            }
            content1.append("%29");
            content1.append("+AND+");
            content1.append(keyValue);
            content1.append("[tiab]");
        } else if (diseaseNameNum == 1) {
            String diseaseName = diseaseNames.get(0);
            //only search this gene
            if (diseaseName.equals("ANY")) {
                content1.append(keyValue);
                content1.append("[tiab]");
            } else if (diseaseName.indexOf(' ') >= 0) {
                //()
                content1.append("%28");
                content1.append(diseaseName.replace(' ', '+'));
                content1.append("[tiab]%29");
                content1.append("+AND+");
                content1.append(keyValue);
                content1.append("[tiab]");
            } else {
                content1.append(diseaseName);
                content1.append("[tiab]");
                content1.append("+AND+");
                content1.append(keyValue);
                content1.append("[tiab]");
            }
        } else {
            return null;
        }


        content1.append("&");
        content1.append(URLEncoder.encode("datetype", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("edat", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("retmax", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("100", "UTF-8"));

        try {
            // URL of CGI-Bin script.
            URL url = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?");

            //http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%28Spinocerebellar+ataxia[mesh]%29+AND+TGM6[tiab]&datetype=edat&retmax=100
            TimeUnit.SECONDS.sleep(4);

            // URL connection channel.
            HttpsURLConnection urlConn = (HttpsURLConnection)url.openConnection();
            // Let the run-time system (RTS) know that we want input.
            urlConn.setDoInput(true);
            // Let the RTS know that we want to do output.
            urlConn.setDoOutput(true);
            // No caching, we want the real thing.
            //urlConn.setUseCaches(false);
            //set timeout time
            urlConn.setConnectTimeout(700000);

            //send data
            // Send quest output.
            OutputStreamWriter wr = new OutputStreamWriter(urlConn.getOutputStream());
            //System.out.println(content1);
            wr.write(content1.toString());
            wr.flush();
            wr.close();

            /*
            BufferedReader rd = new BufferedReader(new InputStreamReader(urlConn.getInputStream()));
            while ((line = rd.readLine()) != null) {
            line = line.trim();
            if (line.length() == 0) {
            continue;
            }
            System.out.println(line);
            }
            rd.close();
             * 
             */


            DocumentBuilder dombuilder = domfac.newDocumentBuilder();
            Document doc = dombuilder.parse(urlConn.getInputStream());
            Element root = doc.getDocumentElement();

            NodeList nodes = root.getElementsByTagName("Count");
            if (nodes.getLength() == 1) {
                Element e = (Element) nodes.item(0);
                Text t = (Text) e.getFirstChild();
                if (t == null) {
                    return "";
                } else {
                    int count = Util.parseInt(t.getNodeValue());
                    if (count == 0) {
                        return "";
                    }
                }
            }
            /*
            nodes = root.getElementsByTagName("IdList");
            if (nodes.getLength() == 1) {
            Node item = nodes.item(0);
            for (Node node = item.getFirstChild(); node != null; node = node.getNextSibling()) {
            if (node.getNodeType() == Node.ELEMENT_NODE) {
            if (node.getNodeName().equals("Id")) {
            String name = node.getFirstChild().getNodeValue();
            System.out.println(name);
            }
            }
            }
            }
             * 
             */


            nodes = root.getElementsByTagName("Id");
            for (int i = 0; i < nodes.getLength(); i++) {
                Element subNode = (Element) nodes.item(i);
                Text t = (Text) subNode.getFirstChild();
                idList.append(t.getNodeValue());
                idList.append(';');
            }

            if (idList.length() > 0) {
                return idList.substring(0, idList.length() - 1);
            } else {
                return "";
            }
        } catch (SocketException ex) {
            System.out.print("... failed!");
            return null;
        } catch (IOException ex) {
            System.out.print("... failed!");
            return null;
        }

    }

    public String dbSNPIDESearch(String[] rsList) throws Exception {
        String line = null;

        int index1 = 0;
        int index2 = 0;

        StringBuilder content1 = new StringBuilder();
        content1.append(URLEncoder.encode("db", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("snp", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("retmode", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("text", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("rettype", "UTF-8"));
        content1.append("=");
        content1.append(URLEncoder.encode("FLT", "UTF-8"));
        content1.append("&");
        content1.append(URLEncoder.encode("id", "UTF-8"));
        content1.append("=");

        StringBuilder idList = new StringBuilder();


        try {
            // URL of CGI-Bin script.
            URL url = new URL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?");

            //http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=3000860&retmode=text&rettype=FLT

            StringBuilder conten = new StringBuilder();

            for (String snp : rsList) {
                // URL connection channel.
                URLConnection urlConn = url.openConnection();
                // Let the run-time system (RTS) know that we want input.
                urlConn.setDoInput(true);
                // Let the RTS know that we want to do output.
                urlConn.setDoOutput(true);
                // No caching, we want the real thing.
                //urlConn.setUseCaches(false);
                //set timeout time
                urlConn.setConnectTimeout(700000);

                //send data
                // Send quest output.
                OutputStreamWriter wr = new OutputStreamWriter(urlConn.getOutputStream());
                //System.out.println(content1);
                wr.write(content1.toString() + URLEncoder.encode(snp.substring(2), "UTF-8"));
                wr.flush();

                BufferedReader rd = new BufferedReader(new InputStreamReader(urlConn.getInputStream()));
                while ((line = rd.readLine()) != null) {
                    line = line.trim();
                    if (line.length() == 0) {
                        continue;
                    }
                    // System.out.println(line);
                    conten.append(line);
                }
                rd.close();
                wr.close();
                if (conten.indexOf("stop-gained") >= 0) {
                    System.out.println(snp + "\tstop-gained");
                } else if (conten.indexOf("stop") >= 0) {
                    System.out.println(snp + "\tstop-lossed");
                } else if (conten.indexOf("missense") >= 0) {
                    System.out.println(snp + "\tmissense");
                } else if (conten.indexOf("splice") >= 0) {
                    System.out.println(snp + "\tsplice");
                } else if (conten.indexOf("synonymous") >= 0) {
                    System.out.println(snp + "\tsynonymous");
                } else {
                    System.out.println(snp + "\tfailed");
                }
                conten.delete(0, conten.length());
                TimeUnit.SECONDS.sleep(3);
            }

            if (idList.length() > 0) {
                return idList.substring(0, idList.length() - 1);
            } else {
                return "";
            }
        } catch (SocketException ex) {
            System.out.print("... failed!");
            return null;
        } catch (IOException ex) {
            System.out.print("... failed!");
            return null;
        }

    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            List<String> dl = new ArrayList<String>();
            /*
            dl.add("Amyotrophic lateral sclerosis");
            dl.add("motor neurone disease");
            dl.add("Lou Gehrig's disease");
             * 
             */
            dl.add("schizophrenia");

            List<String> filter = new ArrayList<String>();
            filter.add("gene");
            filter.add("genes");
            filter.add("mRNA");
            filter.add("protein");
            filter.add("proteins");
            filter.add("transcription");
            filter.add("transcript");
            filter.add("transcripts");
            filter.add("expressed");
            filter.add("expression");
            filter.add("expressions");
            filter.add("locus");
            filter.add("loci");
            filter.add("SNP");
            NCBIRetriever ncbi = new NCBIRetriever();
            /*
            List<String[]> genes = new ArrayList<String[]>();
            int[] indexes = new int[]{1, 4};
            
            // String[] selectedGenes = new String[]{"CALHM1", "MCHR1", "THRB", "TACR2", "BAG4", "ANK3", "CALN1", "DLG1", "MTUS2", "GNL3", "ANAPC1", "NOL4", "STK24", "MAN1A1", "MBD1", "PTPRO", "ECT2", "BICD1", "PRKCB", "XPR1", "RYR2", "PPM1L", "PPM1M", "SNX17", "BTN2A1", "ASAP1", "TMEM219", "ADRBK2", "ZNRD1", "ZNF513", "SORCS3", "PLAA", "RNF165", "CEP170", "BCL11B", "PBRM1", "PRSS16", "GTF3C2", "ZNF804A", "HCN1", "NGEF", "RRP12", "SNX25", "SLC25A12", "FREM3", "ZBED4", "DDX50", "ATP6V0A1", "MIF4GD", "GGA3", "DNM2"};
            // String[] selectedGenes = new String[]{"CYP2D6", "MPV17", "RANGAP1", "ITSN1", "WDR82", "WDR73", "CUL9", "DHX34", "YEATS2", "RPL12", "PTBP2", "LSM1", "DDAH2", "ZNF44", "TTC7B", "KIF5C", "SPIRE2", "PROX2", "TMEM132E", "PTPRO", "ST13", "SERPINF1", "WDR88", "SIPA1L1", "MAPK3", "RYR2", "RPS11", "EIF2AK2", "TP53TG1", "NAGLU", "PUS3", "PPP2R3A", "ZNF79", "GNAI3", "SNORD48", "C18orf21", "NTAN1", "PRRC2A", "PRRC2B", "ZNF513", "JRK", "SLC29A1", "TAOK2", "LY6G5B", "PPAPDC1B", "TMEM110", "ZSCAN2", "USF1", "SFMBT1", "SNORD52", "FXR1", "PPP2R3C", "LIMA1", "BCKDK", "KCNAB2", "DAXX", "YBX2", "DFNB31", "RRN3", "ZNF680", "ZNF737", "TMCO6", "STX4", "PIK3C2A", "SLC22A7", "ANP32E", "EXOSC1", "CLIC1", "BBS10", "LRRC26", "ZSWIM4", "FAM114A2", "XPR1", "ZNF692", "EIF5A2", "PARVB", "ZNF276", "SKIV2L", "C6orf25", "QTRT1", "SF3B1", "DHX16", "BRD1", "PPP1R18", "LRSAM1", "HAPLN4", "STRADA", "TRIM24", "ZNF668", "PKNOX2", "HARBI1", "STAB1", "WHSC1L1", "DGKZ", "DPYD", "KATNAL2", "DUSP7", "CDH10", "LRP4", "DNM1", "DNM2", "GNA13", "C1orf122", "SLC44A2", "DLK2", "CNOT7", "ST3GAL1", "ATAT1", "SMARCD2", "ANK3", "ATF6B", "TUBG2", "FNDC3A", "MTUS2", "NUDT1", "SPTLC1", "STIM2", "CORO7", "HLA-B", "TECR", "BICD1", "EP300", "SGSM1", "WBP2NL", "SMARCA5", "NFE2L1", "MFAP3", "SETD1A", "NKAIN2", "FHDC1", "TTC8", "MDK", "CXXC1", "TUBB", "ISYNA1", "RNF165", "SNF8", "DOPEY1", "STK19", "UPK2", "HCN1", "UCN", "EPB41", "NIN", "SEC11A", "GRIN1", "SHANK3", "S100B", "DDX50", "RPL37A", "IER2", "CAD", "GABBR2", "CCHCR1", "HMOX2", "ALAS1", "TRIM8", "EXOG", "CDK5RAP3", "INO80E", "CASP2", "IK", "AVL9", "NUP85", "PRPF3", "MAN1A1", "ARL3", "PPM1F", "PPM1G", "RFWD3", "ZFYVE28", "TXNRD2", "MFSD10", "PPM1L", "PPM1M", "NCAN", "HTR7P1", "SNX19", "GALNT3", "STON2", "CTF1", "SNX17", "BTN2A1", "BTN2A2", "ATP5G1", "VARS", "FAM13B", "ESPN", "SNX21", "GTF3C2", "PDCD11", "CNTD1", "CCP110", "C1orf173", "GDPD5", "XPNPEP3", "TRANK1", "DBI", "ZBED4", "MIF4GD", "GGA3", "C2orf82", "CALHM2", "VAPA", "HMGN4", "TRMT1L", "CHRNA5", "CREB3L1", "APOM", "C9orf72", "CHRNA3", "HCG27", "ANAPC1", "F11R", "PARP16", "PGM3", "CATSPER2", "MTF1", "ABT1", "FAM177A1", "NEK4", "BTN3A3", "GIGYF2", "GLG1", "KAZN", "CLCN3", "NEK1", "ASAP1", "ADRBK2", "ZNRD1", "CHEK1", "FCGRT", "SFXN2", "SLC32A1", "PLAA", "SIAE", "ACE", "CEP170", "NMRAL1", "PBRM1", "STX10", "ZBTB48", "C21orf91", "CRIP3", "CFB", "NPL", "HNRNPA0", "C6orf1", "DCP1B", "DCP1A", "NOTCH4", "TCF19", "GOSR1", "ALG12", "C4orf36", "ABCF1", "MCHR1", "ZMAT2", "VPS37A", "ZBBX", "FLT3LG", "FANCL", "BAG4", "BAG6", "DIP2B", "SEMA3G", "ACTR5", "ATP5L", "CALN1", "AGAP1", "FANCA", "DNAJC19", "LRRC37B", "RBBP5", "MBD5", "MBD1", "STX1B", "CHMP1A", "ZDHHC16", "TNFRSF10B", "ANKRD24", "DYNC1I2", "TOP3B", "IFT74", "THAP3", "ATG13", "PPP6R1", "NMB", "PPP6R3", "SUFU", "CHADL", "STAT6", "STT3B", "SUMO2", "BCL11B", "SEPSECS", "FYCO1", "SPATA7", "RRP12", "NOS1", "ILF3", "AMBRA1", "RPL7L1", "PBX1", "NOP56", "CD14", "SEPT3", "THRB", "GABRB1", "MRPS30", "TLR9", "TBC1D12", "ROBO1", "CLP1", "ROBO2", "NALCN", "GNL3", "MAGOH", "FOXN2", "SPA17", "NEBL", "TRIM38", "TRIM36", "PSMA4", "MYNN", "NRBP2", "GPN1", "HAX1", "ANAPC13", "COPZ1", "SNX1", "SRF", "RBX1", "GFM2", "EIF3B", "SCARB2", "KAT8", "PRSS53", "TTC21B", "MRPS21", "MUSTN1", "FOXP1", "ATE1", "TRIM54", "PSMC5", "SPCS1", "KCTD13", "C6orf48", "EIF4E2", "SLC25A17", "OPCML", "TACR2", "BZW2", "PHC3", "POU5F1", "POMT1", "ARHGAP1", "ITIH4", "ITIH3", "PLCB2", "CILP", "ADIPOR2", "MRPS7", "DAPK2", "ECT2", "PRKCB", "ARHGAP30", "TYK2", "ACVR2B", "MANEAL", "NRM", "CUX1", "PRKD3", "C11orf24", "TNFRSF25", "HAT1", "DNAH1", "DDHD2", "TYMP", "NPAS3", "NUDT21", "C6orf136", "UCK1", "SCO2", "ZNF605", "ZNF804A", "NGEF", "ZBTB7A", "NDUFA2", "UBE4A", "L3MBTL2", "NDUFA6", "GPANK1", "CWF19L2", "IGSF9B", "ATXN7L2", "RPRD2", "TDP2", "SP4", "SYNM", "TMEM86B", "SLC5A6", "TSSK6"};
            //String[] selectedGenes = new String[]{"ABCA2", "ANKRD11", "CCDC57", "CDC42BPB", "RNFT2", "LRP1", "INPPL1", "TNPO3", "SETD1B", "MACF1", "CEP170", "LRRC10B", "ENPP1", "NAV2", "SH2B3", "PPP6R3", "GALK2", "LDHB", "ABCA5", "TTN", "ADAM17", "CREB5", "CC2D2B", "HTRA2", "ADAMTS15", "LAMC3", "DICER1", "SEMA3F", "HFM1", "SPATA5L1", "USP48"};
            String[] selectedGenes = new String[]{"NAV2", "CDC42BPB", "CCDC57", "SETD1B", "MACF1", "ANKRD11", "RNFT2", "LRP1", "ABCA2"};
            Set<String> selectedGeneSet = new HashSet<String>(Arrays.asList(selectedGenes));
            LocalFile.retrieveData("resources/HgncGene.txt", genes, indexes, "\t");
            
            String pubmedID = null;
            List<String> geneNames = new ArrayList<String>();
            for (String[] item : genes) {
            if (!selectedGeneSet.contains(item[0])) {
            continue;
            }
            geneNames.clear();
            for (String g : item) {
            if (g == null || g.trim().length() == 0) {
            continue;
            } else {
            if (g.indexOf(",") >= 0) {
            String[] names = g.split(",");
            geneNames.addAll(Arrays.asList(names));
            } else {
            geneNames.add(g);
            }
            }
            }
            while ((pubmedID = ncbi.pubMedIDESearch(dl, geneNames, filter)) == null) {
            System.out.print("reconnecting...");
            }
            System.out.println(item[0] + "\t" + item[1] + "\t" + pubmedID);
            }
             * 
             */
            //annovar
         // String[] snps = new String[]{"rs3000860", "rs373006661", "rs4852974", "rs283413", "rs26821", "rs4535533", "rs3015858", "rs111877208"};
          
          // String[] snps = new String[]{"rs201659620", "rs28736716", "rs386670358", "rs544311105", "rs1022477"};
            //String[] snps = new String[]{"rs543816007", "rs7513079", "rs2994114", "rs1830486", "rs2990550", "rs2419525", "rs74630591", "rs28453011", "rs61772344", "rs141295085", "rs2294228", "rs9441132", "rs12059286", "rs12057989", "rs678781", "rs28712116", "rs376108642", "rs114576842", "rs940571", "rs2288072", "rs2001490", "rs10206850", "rs7653834", "rs4602367", "rs2006764", "rs13119659", "rs11548832", "rs324689", "rs10015804", "rs3798130", "rs2740588", "rs17844512", "rs11167743", "rs2860540", "rs11965538",  "rs77504727", "rs9276571", "rs9276572", "rs6913546", "rs80081867", "rs9476080", "rs9476081", "rs62398997", "rs77436138", "rs62398998", "rs146095374", "rs386714630", "rs2960999", "rs6942733", "rs2911690", "rs11003863", "rs73453180", "rs375261215", "rs77976451", "rs2286168", "rs370963766", "rs112485268", "rs74903687", "rs2114566", "rs61736852", "rs11564170", "rs10506156", "rs2933352", "rs1444220", "rs225014", "rs3743242", "rs9921412", "rs7185272", "rs12708923", "rs4788591", "rs199828137", "rs76327478", "rs534522486", "rs7245978", "rs312185", "rs373855535", "rs11880097", "rs1051457", "rs4815025", "rs2236178", "rs1744769", "rs1139793", "rs71296614", "rs132736", "rs1022478"};
          // String[] snps = new String[]{"rs3000859", "rs3010876", "rs141192612", "rs2419526", "rs78540200", "rs58689366", "rs3925698", "rs28691033", "rs1043478", "rs74341264", "rs2020869", "rs2285175", "rs1778541", "rs4665667", "rs2121743", "rs17857100", "rs2353625", "rs3789149", "rs2006748", "rs2241894", "rs31547", "rs2697531", "rs9276571", "rs4406234", "rs556042271", "rs368918103", "rs7830743", "rs11188637", "rs74815769", "rs76287086", "rs73453183", "rs78788561", "rs12099154", "rs11176584", "rs370602169", "rs4788590", "rs556962104", "rs4646523", "rs11673150", "rs6000172", "rs6000173", "rs2227169"};
          //  String[] snps = new String[]{"rs12023156", "rs78095001", "rs2298847", "rs139300"};
            //VEP
          // String[] snps = new String[] {"rs28548017", "rs3000860", "rs373006661", "rs26821", "rs4535533", "rs3015858", "rs111877208", "rs15895", "rs1022220"};
         //String[] snps = new String[]  {"rs117944955", "rs6671527", "rs201659620", "rs28736716", "rs61804988", "rs12568784", "rs3731608", "rs79635065", "rs386670358", "rs191131837", "rs3130453", "rs544311105", "rs7460625", "rs3741132", "rs1043149", "rs202166921", "rs9973206"};
        // String[] snps = new String[]   {"rs199697037", "rs530993551", "rs11260579", "rs1137003", "rs2494620", "rs61744700", "rs10864628", "rs2274327", "rs2274328", "rs61746465", "rs2274333", "rs543816007", "rs1769774", "rs1063776", "rs1830486", "rs1046331", "rs763821", "rs1052576", "rs10927887", "rs2419525", "rs79308175", "rs631357", "rs192833530", "rs143770846", "rs196432", "rs79851686", "rs2297710", "rs3813803", "rs74467253", "rs3903683", "rs12026290", "rs376790279", "rs10379", "rs145091069", "rs199685128", "rs6425977", "rs179472", "rs3738496", "rs2485652", "rs6659553", "rs10789501", "rs61746559", "rs12047252", "rs212991", "rs2494876", "rs537478490", "rs11206407", "rs1655519", "rs1655518", "rs570218", "rs147744184", "rs190299683", "rs1147990", "rs199845652", "rs3008858", "rs1886686", "rs1884444", "rs1133891", "rs273259", "rs34932081", "rs3820093", "rs76005638", "rs1060622", "rs2815413", "rs1801265", "rs4423040", "rs1277207", "rs7483", "rs1058885", "rs2229165", "rs2256721", "rs10776792", "rs149515458", "rs12059286", "rs12057989", "rs678781", "rs2147326", "rs2762745", "rs2798893", "rs2230061", "rs61746299", "rs181187833", "rs200529004", "rs3818831", "rs2105117", "rs2006940", "rs146128753", "rs6668968", "rs6685323", "rs2228145", "rs138911736", "rs4745", "rs61732805", "rs493446", "rs7513351", "rs6668857", "rs857685", "rs140801792", "rs1319080", "rs6427504", "rs4656926", "rs11576830", "rs753856", "rs368151945", "rs913257", "rs12118933", "rs2236410", "rs10911390", "rs1174657", "rs371872781", "rs1174658", "rs16861394", "rs137935651", "rs10801578", "rs554389553", "rs2820289", "rs2292822", "rs2250377", "rs1065761", "rs61745299", "rs35920428", "rs4951168", "rs1150258", "rs11119315", "rs2235371", "rs149508727", "rs940571", "rs2437150", "rs1033325", "rs1055851", "rs2275687", "rs1126627", "rs201959745", "rs4925663", "rs6695302", "rs11583410", "rs6702693", "rs6657127", "rs386641879", "rs6664332", "rs4478844", "rs79716074", "rs6712141", "rs6710817", "rs12623642", "rs4665855", "rs144812943", "rs6715286", "rs2304678", "rs1919125", "rs1919126", "rs1919127", "rs1919128", "rs28381983", "rs3749147", "rs2098386", "rs3213746", "rs4670800", "rs2278586", "rs3816182", "rs3816183", "rs9309107", "rs6720173", "rs698761", "rs115212523", "rs2293275", "rs1045910", "rs2044693", "rs1063479", "rs4453725", "rs2228203", "rs3796097", "rs6748672", "rs10199088", "rs13010513", "rs1813160", "rs531051730", "rs1063588", "rs6707475", "rs6886", "rs6643", "rs1138484", "rs78639210", "rs536823787", "rs2271038", "rs76878531", "rs1044575", "rs138237414", "rs1132267", "rs3731607", "rs567696052", "rs10208844", "rs3905317", "rs1112542", "rs16831235", "rs4954449", "rs7575451", "rs7592429", "rs1432273", "rs3749006", "rs16853344", "rs492594", "rs17544587", "rs4668246", "rs6736609", "rs10185178", "rs6433711", "rs66677602", "rs922984", "rs922985", "rs6729801", "rs288334", "rs13006529", "rs3769823", "rs2043449", "rs2289025", "rs3762568", "rs386654966", "rs2229571", "rs4674361", "rs1127102", "rs1043160", "rs11689281", "rs147759986", "rs2289912",  "rs58277463", "rs6728493", "rs6732185", "rs60617721", "rs62184587", "rs3815305", "rs7424328", "rs6442827", "rs279552", "rs73126218", "rs696217", "rs2276749", "rs12629133", "rs3731130", "rs2228000", "rs7653834", "rs4602367", "rs370892144", "rs156265", "rs1274958", "rs200197381", "rs11079", "rs12495805", "rs6794", "rs7642448", "rs34788938", "rs9834639", "rs370934093", "rs1466685", "rs2071203", "rs2269432", "rs117179004", "rs138951582", "rs9856575", "rs11177", "rs2289247", "rs146224729", "rs2291498", "rs74714110", "rs73840323", "rs73840324", "rs62247157", "rs200756071", "rs200420264", "rs73840339", "rs72503535", "rs200618735", "rs201780087", "rs62250084", "rs62250087", "rs137915102", "rs111230396", "rs202061186", "rs75085872", "rs112669644", "rs77715040", "rs76757483", "rs79811623", "rs80170761", "rs143244597", "rs113197025", "rs113078821", "rs199898382", "rs113494751", "rs112828508", "rs144106206", "rs111677009", "rs78297221", "rs77394685", "rs75759940", "rs76175438", "rs139448820", "rs76263679", "rs62250107", "rs77322475", "rs73843012", "rs74861398", "rs77649281", "rs79870536", "rs375903439", "rs62250110", "rs75467043", "rs142370851", "rs386662584", "rs199839139", "rs78514077", "rs75388892", "rs7653652", "rs4857302", "rs793440", "rs340151", "rs2272022", "rs9826308", "rs3732813", "rs2306858", "rs2291465", "rs200847413", "rs2070179", "rs34006803", "rs2332285", "rs77786095", "rs9833275", "rs1801019", "rs3772809", "rs140693", "rs9883988", "rs367734078", "rs1061409", "rs3816527", "rs28364680", "rs2068178", "rs2290480", "rs2276805", "rs9836841", "rs2304456", "rs710446", "rs7624750", "rs554846023", "rs199822551", "rs201895754", "rs79358488", "rs75737861", "rs373469782", "rs28542401", "rs28438604", "rs372590482", "rs562508581", "rs2911272", "rs413807", "rs430037", "rs56359992", "rs391928", "rs577495599", "rs75588776", "rs7372250", "rs3103959", "rs556018023", "rs71617321", "rs6776064", "rs71611508", "rs61742251", "rs546982384",  "rs1128427", "rs138916724", "rs6811863", "rs2261167", "rs2229940", "rs117401196", "rs141276381", "rs12233719", "rs2306174", "rs13119659", "rs11548832", "rs324689", "rs12507099", "rs6831595", "rs12646270", "rs11722476", "rs7439869", "rs575588990", "rs1126673", "rs1126671", "rs58374767", "rs2715591", "rs4698803", "rs769242", "rs1048201", "rs11098943", "rs541196756", "rs11935573", "rs1352714", "rs11559290", "rs9410", "rs2228119", "rs2463447", "rs4317244", "rs6846627", "rs10035612", "rs10035653", "rs201031593", "rs7728667", "rs386684404", "rs115086471", "rs2966952", "rs2438652", "rs387906317", "rs2278008", "rs34677", "rs2287939", "rs17521570", "rs1035480", "rs2278329", "rs1801033", "rs1499280", "rs2233213", "rs2548612", "rs552295796", "rs2431352", "rs2909888", "rs168939", "rs2366777", "rs4916684", "rs4916685", "rs16868972", "rs890757", "rs25640", "rs7730228", "rs3798130", "rs479632", "rs6890689", "rs6865472", "rs250430", "rs2530245", "rs555552385", "rs2240696", "rs2240695", "rs9686540", "rs11167600", "rs3733709", "rs3733708", "rs7701755", "rs2240694", "rs11167605", "rs3822346",  "rs4141841",  "rs17844301", "rs17844309", "rs3756332", "rs3756331", "rs3733705", "rs6580012", "rs251353", "rs251354", "rs251355", "rs369636", "rs369639", "rs251362", "rs630162", "rs17844348", "rs246074", "rs191251073", "rs31849", "rs13157538", "rs246707", "rs246703", "rs2950845", "rs2950844", "rs17844498", "rs2740582", "rs2697541", "rs17844646", "rs17844664", "rs13176519", "rs549792075", "rs2907323", "rs702386", "rs2910329", "rs2910005", "rs35592458", "rs371627839", "rs7736541", "rs2240697", "rs7709485", "rs15251", "rs6556276", "rs1650893", "rs1211554", "rs576529864", "rs619483", "rs7166", "rs643232", "rs2076299", "rs28763971", "rs9370867", "rs198845", "rs567017585", "rs198844", "rs7766641", "rs141811762", "rs149494377", "rs9368537", "rs112779382", "rs9380030", "rs9986596", "rs12000", "rs1635", "rs9461446", "rs1679709", "rs568537725", "rs6922302", "rs853678", "rs1361385", "rs1416918", "rs2232434", "rs2232430",  "rs386698422", "rs3734815", "rs141466884", "rs757259", "rs3094124", "rs2233980", "rs1265053", "rs142530210", "rs3130981", "rs707913", "rs386698868", "rs2233952", "rs386698869", "rs9501057", "rs1063646", "rs537216792", "rs2073721", "rs3130932",             "rs2269475", "rs1046089", "rs3132453", "rs10885", "rs2242655", "rs3130617", "rs9267532", "rs9267547", "rs2242653", "rs141454939", "rs569405724", "rs2607015", "rs2753960", "rs2227956", "rs12661281", "rs7887", "rs17421133", "rs3134604",  "rs1135216", "rs1057141",   "rs386699870", "rs80330773", "rs1042140",    "rs3130257", "rs3130100", "rs3116713",  "rs1150781", "rs34109614", "rs2239808", "rs11756091", "rs2273962",  "rs35628750", "rs3749903", "rs2242416",  "rs12192544", "rs6913546", "rs80081867", "rs9476080", "rs9476081", "rs62398997", "rs77436138", "rs62398998", "rs592121", "rs3778005", "rs540471794", "rs34012596", "rs17072442", "rs2228547", "rs1064583",  "rs3749895", "rs3828743", "rs6927567", "rs784133",  "rs7742542", "rs62431286", "rs369890576", "rs73558557", "rs181439920", "rs237025", "rs2342767", "rs9371533", "rs9383921", "rs138165204", "rs2151910", "rs541250318", "rs9397449", "rs540825", "rs25683", "rs1128670", "rs683369", "rs2228006", "rs1805323", "rs2640", "rs848016", "rs2074000", "rs2301721", "rs190612258", "rs2302339", "rs2302340", "rs141567309", "rs10216063", "rs531031859", "rs2537188", "rs10899750", "rs2227983", "rs4245575", "rs546721822", "rs17138089", "rs386714630", "rs7056", "rs9955", "rs1128349", "rs28494095", "rs13246460", "rs543957173", "rs2041049", "rs194525", "rs42664", "rs2301680", "rs854524", "rs6962772", "rs2070215", "rs3823646", "rs2075756", "rs184128444", "rs6973420", "rs292501", "rs10265", "rs6464465", "rs10282312", "rs13438232", "rs6967117", "rs201010806", "rs76490935", "rs544887797", "rs201893359", "rs2303929", "rs576660474", "rs17856647", "rs61729490", "rs9693671", "rs56023905", "rs3989701", "rs3177011", "rs149660269", "rs13260331", "rs6990563", "rs2736269", "rs13278959", "rs188668641", "rs17815120", "rs9721012", "rs545142644", "rs75826165", "rs142161778", "rs571781433", "rs2942194", "rs3736032", "rs891398", "rs2472553", "rs578048840", "rs2272736", "rs444772", "rs446227", "rs414352", "rs3739336", "rs959976", "rs9656982", "rs7012186", "rs116726914", "rs144257411", "rs6599528", "rs1735169", "rs6475273", "rs7021572", "rs1071545", "rs74557595", "rs3808869", "rs147165051", "rs2297879", "rs61758536", "rs7044405", "rs79284865", "rs66626885", "rs10821128", "rs1055710", "rs2061634", "rs3747495", "rs3739741", "rs7855745", "rs7869523", "rs12341266", "rs10817493", "rs4837768", "rs4836822", "rs10818503", "rs386738321", "rs1962091", "rs12375547", "rs497632", "rs13298768", "rs574199828", "rs564754", "rs7997", "rs186236478", "rs2296949", "rs2296957", "rs7867616", "rs2073924", "rs886090", "rs12763", "rs945384", "rs945386", "rs2254143", "rs142877518", "rs1052420", "rs553296913", "rs7896053", "rs147888886", "rs145152071", "rs562959917", "rs625223", "rs2230469", "rs1638630", "rs2484180", "rs200934236", "rs1047991", "rs2808096", "rs2505232", "rs12764004", "rs117610424", "rs870957", "rs12269028", "rs41301609", "rs7094610", "rs8187730", "rs11003863", "rs145960427", "rs2452505", "rs3758490", "rs7076156", "rs62620248", "rs10509305", "rs2255607", "rs201068859", "rs906220", "rs2277257", "rs2252996", "rs2487068", "rs1227065", "rs3747862", "rs11558719", "rs111606598", "rs4294502", "rs6480771", "rs11557865", "rs2275047", "rs12782963", "rs369959332", "rs7006", "rs7831", "rs549237264", "rs805722", "rs11192966", "rs3736946", "rs2227309", "rs3750898", "rs3188055", "rs2475298", "rs11200524", "rs139450310", "rs372742362", "rs3781409", "rs7095325", "rs2804003", "rs4369319", "rs1046175", "rs2255246", "rs2686894", "rs1045288", "rs373900649", "rs7107362", "rs7128044", "rs7116130", "rs7128029", "rs1059091", "rs4963198", "rs2740375", "rs7126805", "rs6682", "rs148793744", "rs77438942", "rs75826443", "rs373360588", "rs33988517", "rs76307106", "rs34053383", "rs76461263", "rs77885750", "rs377056834", "rs35840539", "rs34507884", "rs190754396", "rs3741231", "rs2010722", "rs78745657", "rs201073099", "rs11036912", "rs5006888", "rs393044", "rs2017433", "rs2017434", "rs2047456", "rs7935564", "rs75758940", "rs7397032", "rs10769671", "rs1050239", "rs12806437", "rs10839659", "rs12801394", "rs2568023", "rs78813996", "rs7948121", "rs1056951", "rs4243927", "rs2071461", "rs2896586", "rs2071460", "rs2896585", "rs35411689", "rs2271001", "rs1443548", "rs1443549", "rs9666607",  "rs567974104", "rs33952257", "rs141213903", "rs2306029", "rs2306033", "rs6485702", "rs12286721", "rs73453180", "rs375261215", "rs77976451", "rs2063276", "rs2063277", "rs10897403", "rs1905055", "rs201530882", "rs76688946", "rs544864", "rs73481226", "rs55920775", "rs145869982", "rs2229177", "rs542998", "rs186408789", "rs184318096", "rs377667279", "rs17857132", "rs3741378", "rs604630", "rs7947504", "rs377751628", "rs145062279", "rs2282488", "rs4379869", "rs2512526", "rs10501429", "rs494791", "rs7131178", "rs61995873", "rs3765620", "rs693001", "rs659243", "rs2640785", "rs12575909", "rs7124407", "rs2303436", "rs190646973", "rs1356428", "rs3135507", "rs2229113", "rs643788", "rs100803", "rs643423", "rs1815811", "rs78682406", "rs10790715", "rs3088241", "rs549990", "rs622756", "rs2131535", "rs3751037", "rs2532502", "rs714774", "rs562369888",  "rs61733180", "rs1057077", "rs2024301", "rs1860926", "rs1860967", "rs10219561", "rs7308811", "rs11053646", "rs1141715", "rs619381", "rs3741845", "rs200388520", "rs370217783", "rs10845279", "rs10845280", "rs10845281", "rs12226919", "rs12226920", "rs112485268", "rs74903687", "rs11050243", "rs149097385", "rs2933352", "rs11168230", "rs2228570", "rs3741628", "rs12379", "rs7954976", "rs2608009", "rs6580873", "rs79897879", "rs2013335", "rs17101603", "rs2657881", "rs35493121", "rs7397167", "rs2228224", "rs2228226", "rs774895", "rs17012533", "rs9262", "rs12580153", "rs1048906", "rs11109968", "rs11109969", "rs190109616", "rs7969300", "rs1131476", "rs1051042", "rs77534891", "rs1169305", "rs531391976", "rs550483937", "rs371346772", "rs14259", "rs10773019",  "rs12423651", "rs78905764", "rs75085951", "rs559684811", "rs41292167", "rs3864970", "rs35498994", "rs17849654", "rs17646069", "rs1044856", "rs3764147", "rs7998427", "rs1061472", "rs3742289", "rs199765767", "rs1047740", "rs2391191", "rs3742191", "rs41288616", "rs41288620", "rs527819694", "rs2234632", "rs180694682", "rs9624", "rs1263872", "rs1243469", "rs12437266", "rs140996984", "rs34015250", "rs3021119", "rs375638461", "rs7146310", "rs147041732", "rs2295322", "rs1045644", "rs11157729", "rs2073347", "rs12882191", "rs2236316", "rs386777431", "rs146321499", "rs574048486", "rs7142228", "rs370821477", "rs4902176", "rs35949016", "rs544289725", "rs200630198", "rs1133834",  "rs225014", "rs12896583", "rs10141024", "rs2181170",  "rs199531731", "rs113064986", "rs12595504", "rs374230", "rs376432051", "rs73403546", "rs7181742", "rs7180279", "rs2470911", "rs3100139", "rs2576090", "rs8023906", "rs2452524", "rs16963151", "rs17730281", "rs551225", "rs144285622", "rs530593628", "rs2290344", "rs3809530", "rs3809529", "rs61753957", "rs2729835",  "rs2063690", "rs9806123", "rs11071896", "rs7170185", "rs3743093", "rs2277552", "rs2277598", "rs374914417", "rs1048661", "rs3825942", "rs743580", "rs743581", "rs3743211", "rs147106916", "rs76071148", "rs71534208", "rs2271433", "rs149974689", "rs1866928", "rs369194000", "rs7166441", "rs7171194", "rs182240673", "rs7178698", "rs7496668", "rs11247226", "rs183360", "rs531692764", "rs9806942", "rs11248931", "rs448063", "rs2071950", "rs1004041", "rs11248860", "rs1040499", "rs2745103", "rs369792650", "rs1231122", "rs1231123", "rs1139652", "rs2037912", "rs2075639", "rs376032057", "rs2302607", "rs141159056", "rs569963928", "rs4280262", "rs11546303", "rs147546015", "rs76732059", "rs367885695", "rs16967494", "rs7191155", "rs138012549", "rs1140239", "rs2230433", "rs9938550", "rs35713203", "rs11150606", "rs3850114", "rs369576990", "rs8052394", "rs2290182", "rs201330281", "rs1134760", "rs237831", "rs3803650", "rs1131341", "rs74613280", "rs16973716", "rs4788682", "rs567053370", "rs308925", "rs2288020", "rs2288022", "rs2288023", "rs571312243", "rs11557187", "rs2303771", "rs8059048", "rs201164663", "rs28417933", "rs28415940", "rs9921361", "rs55974122", "rs6500437", "rs141278771", "rs2228479", "rs8166", "rs2078478", "rs2272011", "rs2281726", "rs116068630", "rs224496", "rs161400", "rs2976230", "rs1716", "rs2277680", "rs369258303", "rs2286672", "rs188661284", "rs12761", "rs4796555", "rs4562", "rs13290", "rs76686847", "rs3760422", "rs4491591", "rs3803800", "rs1042522", "rs57985356", "rs8522", "rs12453250", "rs9895916", "rs3826543", "rs871841", "rs406446", "rs2074890", "rs200664810", "rs370794931", "rs12603700", "rs1563632", "rs12449313", "rs9912644", "rs77288131", "rs12604020", "rs6561", "rs11567842", "rs4795472", "rs11867457", "rs3760454", "rs887230", "rs58529418", "rs931196", "rs3744355", "rs544654228", "rs7215209", "rs2251660", "rs75117355", "rs142596676", "rs4261597", "rs8074887", "rs9897031", "rs564941164", "rs35612698", "rs1046404", "rs1046403", "rs665268", "rs588703", "rs7218599", "rs5036", "rs149644876", "rs1969205", "rs62077266", "rs62077268", "rs62077270", "rs144398588", "rs79247310", "rs7406910", "rs1133818", "rs3744093", "rs115770790",  "rs1292053", "rs569318874", "rs2727288", "rs2584625", "rs6504459", "rs883541", "rs745143", "rs904384", "rs904383", "rs745142", "rs1566286", "rs35682745", "rs35999669", "rs3744215", "rs6501741", "rs200798799", "rs3744032", "rs3803739", "rs77585764", "rs139139850", "rs2071214", "rs8074547", "rs2066964", "rs34367357", "rs147954634", "rs8084295", "rs1965665",  "rs11081510", "rs2240906", "rs662515", "rs376164368", "rs1049683", "rs548027190", "rs370143194",  "rs563079733", "rs2510019", "rs1942418", "rs2279269", "rs572020",  "rs3911730", "rs3764645", "rs12459404", "rs10417628", "rs735911", "rs352493", "rs61739549", "rs966384", "rs17851959", "rs17851960", "rs475002", "rs2547064", "rs1867691", "rs2230752", "rs5498", "rs1056538", "rs2228615", "rs563676718", "rs7258368", "rs11669628", "rs78161395", "rs609636", "rs609290", "rs11666267", "rs1059519", "rs1059369", "rs1007160", "rs12150889", "rs4805825", "rs112886802", "rs3816032", "rs2290649", "rs1688016", "rs1688005", "rs12110", "rs909072", "rs4806163", "rs12460932", "rs34562867", "rs376496955", "rs382789", "rs3745775", "rs4802029", "rs28512414", "rs17856467", "rs56854837", "rs11670988", "rs10425169", "rs7249222", "rs10853751", "rs284662", "rs1043413", "rs75257969", "rs7245978", "rs11883278", "rs143377407", "rs2302421", "rs454301", "rs439676", "rs435590", "rs365745", "rs366111", "rs4239529", "rs375784428", "rs3746319", "rs372853759", "rs12151338", "rs2571174", "rs203710", "rs35385129", "rs13306210", "rs369933487", "rs5167", "rs735482", "rs762562", "rs2336219", "rs3212986", "rs147675790", "rs312185", "rs373855535", "rs2043211", "rs2302951", "rs2302947", "rs66877153",  "rs8110220", "rs149205469", "rs3786734", "rs374140759", "rs186208217", "rs2303757", "rs919364", "rs2303759", "rs13394", "rs111846329", "rs8109661", "rs3745486", "rs201377040", "rs113934432", "rs386810263", "rs55735528", "rs4801853", "rs4802741", "rs111846937", "rs73049612", "rs10409962", "rs2305772", "rs61736479", "rs75421648", "rs76466205", "rs75090024", "rs8107444", "rs57548937", "rs13382164", "rs2708743", "rs36657", "rs383369", "rs10423424", "rs618835", "rs1051457", "rs605219", "rs1142881", "rs533930457", "rs643347", "rs1049150", "rs45551936", "rs1130513", "rs201009466", "rs1054940",  "rs80009430", "rs3760849", "rs10409531", "rs34051133", "rs1027392", "rs3795003", "rs11879465", "rs144181716", "rs9575", "rs1135196", "rs72484087", "rs142972068", "rs200055053", "rs541896402", "rs6051545", "rs1545", "rs1547", "rs145045986", "rs1205193", "rs536798351", "rs2255255", "rs2255258", "rs2273056", "rs2273058", "rs6138482", "rs2076559", "rs910397", "rs2425049", "rs3752290", "rs2235592", "rs190116287", "rs13037900", "rs6032538", "rs4638862", "rs2275801", "rs397689004", "rs6025606",  "rs6070697", "rs944895", "rs6089219", "rs187680217", "rs114867526", "rs6090457", "rs2226548", "rs2298437", "rs9305426", "rs2507733", "rs574981993", "rs6517549", "rs2297270", "rs220128", "rs220130", "rs2248490", "rs10595", "rs2020945", "rs233239", "rs233252", "rs377549", "rs452472", "rs2020221", "rs944419", "rs363877", "rs446817", "rs369720", "rs8127342", "rs9980129", "rs4818952", "rs587623883", "rs9306111", "rs35548026", "rs17183290", "rs572275194", "rs8190315", "rs7575", "rs4680", "rs165815", "rs426938", "rs530867184", "rs377384823", "rs3177243", "rs62231981", "rs200145797", "rs1894704", "rs1894706", "rs5752330", "rs2014410", "rs713998", "rs538054152", "rs9625679", "rs743920", "rs542860875", "rs151177194", "rs2015035", "rs2072193", "rs11107", "rs1053593", "rs132653", "rs71296614", "rs132736", "rs2075939", "rs855791", "rs9610775", "rs2285178", "rs2076371", "rs137055", "rs13815", "rs368506169", "rs2272837", "rs140523", "rs3764880", "rs12688514", "rs151281801", "rs1047037", "rs1047034", "rs6624595", "rs1045686", "rs2642219", "rs5951328", "rs1298577", "rs147613758", "rs41301507", "rs176024", "rs176026", "rs1139916", "rs2266890", "rs1059703", "rs12012747"};
       //   String[] snps = new String[]   {"rs2307130", "rs78095001", "rs13394619", "rs17850855", "rs2229638", "rs3735887", "rs549545181", "rs10899795", "rs2521998", "rs1671926", "rs10774671", "rs3803354", "rs2958059", "rs2298847", "rs8107859", "rs2074071", "rs139300"};
      //   String[] snps = new String[]    {"rs4970441", "rs2298213", "rs17568", "rs10907179", "rs307349", "rs819980", "rs7533", "rs10797440", "rs61734742", "rs61735051", "rs2274330", "rs370819293", "rs11121691", "rs2275527", "rs2066470", "rs145408535", "rs193921033", "rs3000859", "rs3010876", "rs1737157", "rs1063779", "rs141192612", "rs355025", "rs2294490", "rs2236215", "rs2863452", "rs2419526", "rs2293920", "rs1737428", "rs6692016", "rs4704", "rs12046928", "rs1803531", "rs2294520", "rs386629730", "rs200760163", "rs2273270", "rs3813804", "rs75060926", "rs1050663", "rs1061770", "rs2228550", "rs744455", "rs2270537", "rs3737995", "rs6425816", "rs911218", "rs768000", "rs2242427", "rs7065", "rs6524", "rs41267037", "rs2429050", "rs201492855", "rs2298005", "rs2275276", "rs1250", "rs374946861", "rs1048771", "rs12161316", "rs15921", "rs4926650", "rs9332417", "rs2304313", "rs2304312", "rs3765018", "rs3765017", "rs792310", "rs1051122", "rs12703", "rs786906", "rs2783499", "rs6583048", "rs143563718", "rs2784140", "rs2229155", "rs17030651", "rs2306940", "rs200652089", "rs28691033", "rs1043478", "rs71664005", "rs2762867",  "rs74443150", "rs10749658", "rs17850000", "rs188568243", "rs7172", "rs1196455", "rs12022217", "rs6679449", "rs36107483", "rs1760795", "rs375522819", "rs11264295", "rs2242195", "rs6680600", "rs1142287", "rs1052176", "rs914616", "rs6427322", "rs6668178", "rs386635740", "rs3827757", "rs2427808",  "rs34749866", "rs143083729", "rs33941127", "rs72633678", "rs404508", "rs74341264", "rs17852003", "rs41269662", "rs7531125", "rs2020869", "rs2285175", "rs3738423", "rs11332", "rs6677840", "rs10911392", "rs12045762", "rs2296713", "rs400344", "rs3729547", "rs2228079", "rs59861987", "rs17851127", "rs4950979", "rs569896475", "rs167082", "rs3736963",  "rs1141714", "rs1136448",  "rs61729123", "rs4149229", "rs4837", "rs1051038", "rs3738402", "rs138886515", "rs17846600", "rs10924790", "rs2275685", "rs1885532", "rs4518892", "rs1341864", "rs1341863", "rs10802990", "rs10927387", "rs12028142", "rs6702695", "rs532917909", "rs10170348", "rs16867251", "rs2302942", "rs34080125", "rs7593864", "rs7423300", "rs1402962", "rs2289358", "rs1275536", "rs3811644", "rs2710623", "rs546203668", "rs12712859", "rs149597004", "rs566207204", "rs805412", "rs805423", "rs6740641", "rs848291", "rs12479056", "rs14170", "rs778155", "rs150595565", "rs329498", "rs329497", "rs8827", "rs1137930", "rs2287328", "rs2303606", "rs7210", "rs61735750", "rs34617503", "rs1078004", "rs1009", "rs375756032", "rs11891495", "rs17027011", "rs375321327", "rs5008284", "rs202061207", "rs1519654", "rs11538197", "rs201819489", "rs4848197", "rs3761702", "rs3748914", "rs2276561", "rs1210133", "rs147906017", "rs1126839", "rs10864986", "rs2592595",  "rs4662775", "rs4662594", "rs2290111", "rs1699", "rs1052319", "rs7340193", "rs6423208", "rs11888766", "rs11899779", "rs10928526", "rs6719488", "rs2236783", "rs1525579", "rs4664114", "rs17642086", "rs4556933", "rs13407368", "rs2302694", "rs13385825", "rs2161916", "rs17497636", "rs6751520", "rs6761297", "rs1878583", "rs79541997", "rs1435573", "rs10803917", "rs922986", "rs1143674", "rs188628120", "rs10180793", "rs288280", "rs192618", "rs397843842", "rs45463799", "rs17854823", "rs788023", "rs8539", "rs1050347", "rs80282360", "rs3731700", "rs849563", "rs59358632", "rs3820900", "rs145501967", "rs2287599", "rs2070093", "rs17501837", "rs4324314", "rs2276631", "rs374891392", "rs4674354", "rs3731881", "rs1127101", "rs2293072", "rs1109866", "rs17847405", "rs1976616", "rs13005918", "rs570989435", "rs2305138", "rs3816334", "rs146073833", "rs2645774", "rs2646260", "rs77397270", "rs13382446", "rs4292120", "rs12614632", "rs144572631", "rs3208142", "rs6771714", "rs1705805", "rs112887807", "rs72492998", "rs72624420", "rs544377883", "rs14080", "rs75326306", "rs61729306", "rs7614776", "rs2012153", "rs1274957", "rs144215660", "rs12487911", "rs339697", "rs147924427", "rs4683310", "rs13062723", "rs4858798", "rs9713651", "rs3020779", "rs2291542", "rs1061474", "rs1138536", "rs6769055", "rs71309963", "rs3732506", "rs17058639", "rs186956439", "rs4264746", "rs10489", "rs61188513", "rs6549590", "rs148233498", "rs62247159", "rs62247960", "rs200618735", "rs111381670", "rs372539402", "rs79902440", "rs141393015", "rs76396775", "rs139448820", "rs62250105", "rs79560192", "rs73843011", "rs62250109", "rs145802847", "rs77688099", "rs141238251", "rs79540623", "rs2279720", "rs17852198", "rs11353", "rs1163439", "rs62276974", "rs72491121", "rs369042607", "rs1131265", "rs368245262", "rs150627065", "rs3772126", "rs2289843", "rs56047676", "rs2291078", "rs13146", "rs7613944", "rs58595496", "rs140414363", "rs2659690", "rs2981026", "rs1680778", "rs1979529", "rs876755", "rs11920311", "rs322113", "rs11549806", "rs6794496", "rs2170309", "rs1802904", "rs6806535", "rs150765088", "rs1463725", "rs201104179", "rs10936599",  "rs7622479", "rs7636910", "rs2228291", "rs956732", "rs13069661", "rs7652597", "rs1047115", "rs11538612", "rs187868", "rs9821880", "rs1675943", "rs2432527", "rs74201433", "rs376353447", "rs377407201", "rs74957513", "rs76887830", "rs71635078", "rs71635079", "rs551802747", "rs2948677", "rs74422441", "rs3103955", "rs369770584",  "rs3103958", "rs79482340", "rs77382269", "rs79661483", "rs17852687", "rs61179255", "rs12636891", "rs1147240", "rs573708", "rs1056664", "rs60457064", "rs17857100", "rs2353625", "rs3789149", "rs2006748", "rs4045481", "rs2305185", "rs2073504", "rs61746000", "rs573685829", "rs3765121", "rs567085813", "rs2125313", "rs2271173", "rs1397548", "rs28365063", "rs13039", "rs4859572", "rs11097244", "rs11558468", "rs11558469", "rs12649750", "rs710834", "rs2306802", "rs6823404", "rs1126670", "rs2032349", "rs199992363", "rs7689008",  "rs2276959", "rs34561493", "rs3756122", "rs148074622", "rs62331900", "rs6848883", "rs12504074", "rs34405980", "rs12502935", "rs1131553", "rs193196784", "rs4635850", "rs75112782", "rs6555055", "rs6865765", "rs4975629", "rs7704058", "rs200341614", "rs7724858", "rs162030", "rs367602123", "rs2578617", "rs61748198", "rs6451206", "rs267766", "rs1531545", "rs2233216", "rs1875506", "rs1051846", "rs158921", "rs182190", "rs469039", "rs200848258", "rs4532349", "rs308365", "rs188169797", "rs950692", "rs293035", "rs9667", "rs7724759", "rs62370437", "rs11748794", "rs10078748", "rs376004978", "rs7522", "rs2228112", "rs254286", "rs4835678", "rs10073922", "rs10117", "rs201846062", "rs1800882", "rs2251860", "rs3733707", "rs3822349", "rs3822347", "rs11958868", "rs7710794", "rs251356", "rs251363", "rs251372", "rs251377", "rs251380", "rs155361", "rs155820", "rs17844387", "rs13161592", "rs17844397", "rs17844551", "rs597064", "rs610836", "rs622424", "rs114149224", "rs2910330", "rs2910331", "rs17097300", "rs2075661", "rs1423148", "rs11167756", "rs138061122", "rs2434322", "rs258819", "rs3764930", "rs1042720", "rs57292511", "rs4841", "rs1048723", "rs77739145", "rs192192", "rs4704870",  "rs11546512", "rs2305727", "rs2305728", "rs371742098",  "rs4868663", "rs79403503",  "rs4631", "rs1130857", "rs335435", "rs335438", "rs11541392", "rs4935", "rs148278350", "rs6877400", "rs139672035", "rs375095348", "rs1016835", "rs397516967", "rs549912670", "rs2277105", "rs126405", "rs3757261", "rs9467583", "rs9358871", "rs9358890", "rs2071299", "rs3752419", "rs35549700", "rs1059486", "rs3734541", "rs3800302", "rs3800303", "rs7756481", "rs200981", "rs16868094", "rs200948", "rs200973", "rs904142", "rs201692836", "rs2859348", "rs551240433", "rs200993279", "rs2074482", "rs200942796", "rs2057727", "rs2074474", "rs2074505", "rs1136719", "rs1419693", "rs2286656", "rs1265055", "rs1062470", "rs3094215", "rs3130983", "rs130077", "rs2073722", "rs9501063", "rs2308604", "rs386699190",  "rs386699191", "rs1129640", "rs11229", "rs9366785", "rs14365", "rs538733367", "rs707926", "rs535586",  "rs117905900",  "rs186853831", "rs3208181", "rs2071471",        "rs9277965", "rs3119025", "rs45560246", "rs386700146", "rs2071809", "rs2296934", "rs2274587", "rs2284922", "rs11753141", "rs9394831", "rs201760202", "rs3749904", "rs4711738", "rs61733494", "rs9472022", "rs941849", "rs13296", "rs484757", "rs516582", "rs325008", "rs498512", "rs12528232", "rs3757241", "rs9296573", "rs2038149", "rs4406234", "rs542948", "rs6917467", "rs61757811", "rs62431288", "rs2942", "rs9373491", "rs3924871", "rs3798761", "rs3798763", "rs17079029", "rs2747662", "rs675026", "rs562859", "rs1032141", "rs1571766", "rs927718", "rs1555774", "rs3464", "rs1803989", "rs1867351",  "rs9460113", "rs1042327", "rs1056568", "rs59438885", "rs376896309", "rs1805319", "rs12532895", "rs2639", "rs2230263", "rs929509", "rs3735231", "rs15775", "rs7779633", "rs61756142", "rs2074566", "rs2158218", "rs2301720", "rs386711774", "rs1584614", "rs17170223", "rs10265207", "rs2595701", "rs8580", "rs2293106", "rs75253453", "rs10951936", "rs12669559", "rs2230197", "rs6460052", "rs8565", "rs187692176", "rs7953", "rs74934684", "rs3747807", "rs6972561", "rs7797834", "rs12536287", "rs2301629", "rs12728", "rs114658645", "rs883403", "rs1859690", "rs3823641", "rs12267", "rs3800951", "rs3735241", "rs1043915", "rs863450", "rs1054391", "rs13241786", "rs73163792", "rs73163793", "rs67377634", "rs4520097", "rs1043615", "rs3735642", "rs1053124", "rs397954408", "rs8043", "rs2305324", "rs13232207", "rs1722883", "rs2696877", "rs35360282", "rs6967301", "rs10255702", "rs1804527", "rs4595035", "rs17171113", "rs917124", "rs3800782", "rs2566514", "rs2608293", "rs2487154", "rs60299355", "rs59912051", "rs370818556", "rs1063523", "rs2916747", "rs201656196", "rs3020221", "rs6559167", "rs76373071", "rs2572403", "rs8417", "rs11555738", "rs2285291", "rs4705", "rs2306825", "rs2280444", "rs76263557", "rs200092345", "rs1126677", "rs189318060", "rs2294128", "rs1010156", "rs2294127", "rs2294133", "rs10992", "rs6996616", "rs1812594", "rs7005936", "rs1879182", "rs1879181", "rs2565061", "rs4149253", "rs240951", "rs4733961", "rs139284078", "rs1058720", "rs441800", "rs3739334", "rs3739337", "rs370572629", "rs7814198", "rs1063045", "rs1129152",  "rs11552577", "rs5020517", "rs559932736", "rs4870887", "rs4483140", "rs6470219", "rs6470220", "rs3812471", "rs72713066", "rs12547243", "rs2069568", "rs441914", "rs2304284", "rs2304283", "rs1045248", "rs6558341", "rs2290416", "rs896954", "rs147489825", "rs4875060", "rs4977191", "rs11136256", "rs4317614", "rs2721172", "rs9004", "rs2785333", "rs106033", "rs33962342", "rs200486111", "rs7875404", "rs2228416", "rs3739548",  "rs3739681", "rs11790577", "rs2275422",  "rs7158", "rs73451725", "rs2809270", "rs4145894", "rs7859201", "rs11144089", "rs2228173", "rs45515192", "rs12683119", "rs3747504", "rs3747496", "rs577984845", "rs11581", "rs7472", "rs150839166", "rs1134931", "rs3739696", "rs3810906", "rs12350531", "rs1826232", "rs77100552", "rs2501727", "rs3750494", "rs687434", "rs466994", "rs1854706", "rs2297866",  "rs522328", "rs570570464", "rs482095", "rs11545664", "rs6781", "rs947624", "rs116991188", "rs4837291", "rs2280844", "rs2287363", "rs117012336", "rs10901065", "rs3739494", "rs7862221", "rs2285486", "rs1055432", "rs77905", "rs3827848", "rs2228560", "rs12684650", "rs3812547", "rs59902911", "rs10781510", "rs3812580", "rs11849", "rs9987876", "rs4880091", "rs28407527", "rs575311659", "rs1052333", "rs1129614", "rs10466280", "rs6686", "rs12415754", "rs147288862", "rs11015742", "rs17756919", "rs11597888", "rs74509719", "rs386743529", "rs397730068", "rs10733838", "rs1162753", "rs10509306", "rs34040486", "rs1084004", "rs3747866", "rs1248634", "rs3814213", "rs4244947", "rs376853320", "rs2234978", "rs3758499", "rs200119542", "rs55843714", "rs726176", "rs727852", "rs3750715", "rs2863095", "rs2282295", "rs369670140", "rs2254537", "rs117905984", "rs10749291", "rs2289306", "rs72839700", "rs34265990", "rs2247705", "rs2474328", "rs2767434", "rs2995326", "rs3008334", "rs1046178", "rs2294081", "rs12252", "rs1134578", "rs75782220", "rs12628", "rs7930569", "rs80335370", "rs10902224", "rs151061093", "rs34281988", "rs548901689", "rs74202058", "rs61869004", "rs374545453", "rs3741232", "rs800342", "rs2010719", "rs77135024", "rs10769023", "rs872751", "rs3751006", "rs199733222", "rs4486660", "rs2288283", "rs1043388", "rs1043390", "rs2292195", "rs1128396", "rs571404172", "rs17847539", "rs10839976", "rs6735", "rs2568068", "rs201862125", "rs1056963", "rs376178332", "rs12794714", "rs1140047", "rs11605097", "rs57253324", "rs1442710", "rs6483642", "rs2953310", "rs2434478", "rs20542", "rs2290883", "rs74815769", "rs76287086", "rs73453183", "rs78788561", "rs948847", "rs950803", "rs950802", "rs139600243", "rs2905692", "rs1800007", "rs1109748", "rs1800009", "rs374829632", "rs4693", "rs10897461", "rs614397", "rs947939", "rs4379854", "rs138590813", "rs612448",  "rs1546532", "rs2236684", "rs2004649", "rs151110979", "rs637571", "rs141538682", "rs1784030", "rs11110", "rs545009", "rs486584", "rs200338165", "rs373309307", "rs2471829", "rs9344", "rs1061328", "rs17853270", "rs4944052", "rs1125609", "rs4479014", "rs10751296", "rs2658797", "rs593690", "rs1150360", "rs75746929", "rs61995874", "rs471933", "rs377224881", "rs186237230", "rs470558", "rs1052313", "rs555367", "rs144543614", "rs1176713", "rs2256111", "rs7104819", "rs512703", "rs2511841", "rs587985", "rs1893261", "rs2156634", "rs12099154", "rs373599590", "rs76579229", "rs3751036", "rs4766334", "rs2907608", "rs17845948", "rs216902", "rs740850", "rs1043262", "rs1065691", "rs1803621", "rs11575110", "rs1051409", "rs61733179", "rs1860927", "rs3764021", "rs11054142", "rs11054143", "rs1805522", "rs61754407", "rs369111462", "rs199732157", "rs4963793", "rs114192162", "rs34624361",  "rs2269828", "rs836180", "rs1129406", "rs1049467", "rs697634", "rs4761856", "rs4761784", "rs401926", "rs12818575", "rs7135148", "rs200634786", "rs8916", "rs1546808", "rs4759315", "rs1249378", "rs1136082", "rs2271189", "rs1800194", "rs1800139", "rs1800141", "rs1800154", "rs7308698", "rs1140648", "rs10783816", "rs2229717", "rs11537654", "rs923829", "rs11180483", "rs366527", "rs1515565", "rs3794261", "rs2303626", "rs6615", "rs146235373", "rs1063936", "rs104895337", "rs10849900", "rs10774648", "rs17849837", "rs7174", "rs76543640", "rs2258227", "rs7965804", "rs10847910", "rs7966867", "rs7135542", "rs145795080", "rs2293512", "rs74727297", "rs201788230", "rs3864971", "rs61947037", "rs7981616", "rs111294919", "rs61729909", "rs3764108", "rs4943266", "rs2772364", "rs1127446", "rs2231332", "rs3812896", "rs3014960", "rs9316179", "rs7337140", "rs6313", "rs61749884", "rs3742291", "rs9565152", "rs138253916", "rs12261", "rs2230343", "rs2230342", "rs9555726", "rs1290177", "rs4150731", "rs8809", "rs2275007", "rs201385630", "rs3762158", "rs1139130", "rs1957374", "rs1061040", "rs2295682", "rs178640", "rs2069542", "rs2229661", "rs2273913", "rs1885711", "rs45484197", "rs2025258", "rs4247001", "rs17105087", "rs11558738", "rs2297995", "rs1983764", "rs2073348", "rs2073349", "rs17853629", "rs2273431", "rs3818186", "rs904059", "rs3737171", "rs2296921", "rs7161127", "rs1256049", "rs1741487", "rs2229677", "rs734028", "rs8014577", "rs2075025", "rs2301345", "rs2230237", "rs20578", "rs8017642", "rs3742764", "rs7250", "rs11159286", "rs176960", "rs12434329", "rs2295660", "rs7148456", "rs10135507", "rs1132975", "rs894039", "rs3803322", "rs2816606", "rs1008628", "rs2239669", "rs570914443", "rs1133642", "rs674155", "rs8208", "rs7321", "rs73407109", "rs1197689", "rs1672466", "rs542271554", "rs677895", "rs677845", "rs12594325", "rs7167392", "rs2278857", "rs12438025", "rs3097773", "rs1053492", "rs200642821",  "rs61756712", "rs4774549", "rs6416452", "rs2069133", "rs3809528", "rs332259", "rs3759785", "rs1065080", "rs311889", "rs374524755", "rs12914333", "rs8023358", "rs11634630", "rs1130741", "rs1128933", "rs3743075", "rs8040868", "rs906439", "rs8031107", "rs118149048", "rs7183618", "rs7169981", "rs2301831", "rs6226", "rs17854846", "rs11247253", "rs561563829", "rs11647490", "rs11649031", "rs763053", "rs370024606", "rs237674", "rs2076443",  "rs2235648", "rs2302175", "rs2286466", "rs8460", "rs1046502", "rs117291307", "rs2301802", "rs10641", "rs12032", "rs224206", "rs224207", "rs224208", "rs37811", "rs2270494", "rs1059857", "rs555198215", "rs758044", "rs9673735", "rs1049208", "rs369409228", "rs1049206", "rs1047732", "rs1382390", "rs1970817", "rs201575210", "rs8191328", "rs1050163", "rs1050162", "rs2075511", "rs3743690", "rs1692729", "rs386789922", "rs376468731", "rs11401", "rs61747536", "rs1131543", "rs2450399", "rs11544328", "rs4787643", "rs3751847", "rs9926903", "rs1549295", "rs1549294", "rs1799917", "rs4937", "rs201477611", "rs3743556", "rs2288012", "rs369106170", "rs3607", "rs11540994", "rs372701071", "rs12932850", "rs2303225", "rs1049794", "rs1050361", "rs1050362", "rs10852515", "rs2303275", "rs4887772", "rs2287990", "rs2232503", "rs2232504", "rs577294720", "rs28370522", "rs2288019", "rs2230126", "rs192318511", "rs2306049", "rs8043637", "rs71395334", "rs4785627", "rs28638280", "rs60437616", "rs61741725", "rs62000377", "rs2279349", "rs11649210", "rs1800358", "rs11648433", "rs549670238", "rs61748635", "rs2272012", "rs3809872", "rs8065251", "rs8077638", "rs6828", "rs8067660", "rs887387", "rs2020118", "rs1876444", "rs1050997", "rs3851", "rs1132446", "rs238230", "rs11209", "rs1071705", "rs2301739", "rs1071648", "rs8078571", "rs14309", "rs3744396", "rs11078663", "rs33979567", "rs560435670", "rs1127440", "rs2228129", "rs3803798", "rs4968189", "rs1057086", "rs148339522", "rs149371094", "rs61179727", "rs1242489", "rs3829956", "rs4393623", "rs916823", "rs7216", "rs2242345", "rs1129506", "rs9893935", "rs372804994", "rs58654829", "rs1053651", "rs8073323", "rs199500732", "rs8182306", "rs371064123", "rs1128966", "rs4796712", "rs2070106", "rs2230326", "rs200318623", "rs690941", "rs13229", "rs2070835", "rs3815076", "rs2071167", "rs5910", "rs12164", "rs63750222", "rs62077265", "rs1058201", "rs62077269", "rs11079804", "rs7406798", "rs7207842", "rs2269772", "rs1064055", "rs1057068", "rs34818467", "rs45511291", "rs6503905", "rs11368", "rs968719", "rs13030", "rs2070776", "rs10127", "rs7105", "rs556938635", "rs3897753", "rs1127737", "rs8066909", "rs17853024", "rs2087718", "rs17844857", "rs35857535", "rs9988", "rs4789164", "rs368567805", "rs8669",  "rs62088214", "rs2665998", "rs17845783", "rs142401033", "rs2289527", "rs368028005", "rs56407805", "rs2304854", "rs61757652", "rs7216493", "rs4889848", "rs1139405", "rs376921756", "rs9303026", "rs1056534", "rs1140458", "rs1049684", "rs2278463", "rs8084554", "rs61104666", "rs7228084", "rs2276186", "rs17852712", "rs17852674", "rs319438", "rs8092336", "rs485613", "rs7238987", "rs7230037", "rs41551212", "rs199912442", "rs568504352", "rs77809401", "rs370122679", "rs4806908", "rs4806909", "rs12973948", "rs881768", "rs2074454", "rs556962104", "rs12488", "rs61735591", "rs12984675", "rs546001063", "rs36526", "rs888932", "rs1138253", "rs4807583", "rs10811", "rs549601872", "rs348362", "rs181532763", "rs2252675", "rs3826783", "rs1559186", "rs202180663", "rs2305791", "rs3745255", "rs2304086", "rs12972621", "rs1078264", "rs11085822", "rs10409785", "rs1127307", "rs78992533", "rs2190686", "rs2074262", "rs141533891",  "rs755123", "rs146783459", "rs6743", "rs1804826", "rs964132", "rs1007161", "rs11084673", "rs12150890", "rs10411735", "rs2304102", "rs7258185", "rs3745975", "rs2112800", "rs1056041", "rs3745973", "rs570908694", "rs143821677", "rs2239945", "rs231228", "rs35297478", "rs3761096", "rs61742664", "rs11538454", "rs2301734", "rs4801861", "rs11544093", "rs2230694", "rs2233156", "rs284660", "rs284661", "rs284663", "rs3752172", "rs851609", "rs7255164", "rs7256456", "rs1065178", "rs374110942", "rs453640", "rs7508149", "rs4508518", "rs3746321", "rs3200505", "rs204481", "rs3786505", "rs204468", "rs8105198", "rs13436", "rs35859451", "rs386810061", "rs8108468", "rs2230267", "rs1056917", "rs1056914", "rs6521", "rs3745296", "rs354021", "rs74863643", "rs2305915", "rs283525", "rs16981617", "rs61751956", "rs17849961", "rs73932617", "rs2075802", "rs10413435", "rs9917044", "rs4803129", "rs12460528", "rs366337", "rs7257187", "rs3745420", "rs3813148", "rs78691806", "rs649216", "rs267605676", "rs200114273", "rs652188", "rs79401710", "rs1671151", "rs1077806", "rs977070", "rs61732213", "rs151198103", "rs12978696", "rs3760850", "rs8107758", "rs145011", "rs574516590", "rs3764574", "rs862704", "rs2074070", "rs3746219", "rs11880050", "rs3746218", "rs3746217", "rs10853905", "rs10418774", "rs2305120", "rs3816328", "rs7271033", "rs944110", "rs141972559", "rs45471597", "rs9305125", "rs561231326", "rs563205671", "rs8958", "rs75415902", "rs1178016", "rs8362", "rs9101", "rs16994453", "rs17852625", "rs16991547", "rs660710", "rs28455091", "rs2026022", "rs1205194", "rs538749857", "rs6035051", "rs2273057", "rs138250783", "rs2233832", "rs4911536", "rs200631265", "rs2424922", "rs1998233", "rs6121015",  "rs4812332", "rs16987712", "rs2235611", "rs11274", "rs2076026", "rs13217", "rs707577", "rs6032537", "rs3746493", "rs1537028", "rs73598392", "rs7121", "rs163777", "rs16985970", "rs112350373", "rs2294995", "rs368678392", "rs2273487", "rs2273506", "rs6122130", "rs56130722", "rs11551685", "rs6090457", "rs2738787", "rs201632478", "rs388707", "rs20572", "rs2835655", "rs199919585", "rs220109", "rs4148972", "rs4148973", "rs872331", "rs1133779", "rs9975113", "rs435424", "rs448908", "rs1211103", "rs388827", "rs389167", "rs437149", "rs437334", "rs370092", "rs8126553", "rs8131142", "rs9979457", "rs1051367", "rs200553372", "rs10432965", "rs4819925", "rs5994165", "rs1034859", "rs4488761", "rs455072", "rs4633", "rs4818", "rs769224", "rs17819104", "rs4675", "rs178269", "rs13054014", "rs12484060", "rs2072775", "rs6005977", "rs2074739", "rs11635", "rs3827346", "rs2235321", "rs4820313", "rs147933419", "rs71948", "rs137831", "rs1799932", "rs133369", "rs133383", "rs35318967", "rs138217", "rs2272836", "rs2272838", "rs35820251", "rs1555048", "rs2066770", "rs55861809", "rs1129880", "rs2072876", "rs570272844", "rs743615", "rs128941", "rs12148",  "rs5924530", "rs3213451", "rs2071311", "rs2856733", "rs12687163", "rs6610447", "rs6520683", "rs1023065", "rs4898", "rs11551797", "rs3048", "rs1128363", "rs1264013", "rs2516023", "rs1199470", "rs2296542", "rs2642218", "rs145959014", "rs5974570", "rs74509382", "rs34606958", "rs5983916", "rs16995747", "rs2233055", "rs5201", "rs13397"};
         String[] snps = new String[] {"rs78453011", "rs2419525", "rs12026633", "rs78667364", "rs17692792", "rs9845816", "rs4685171", "rs6950750", "rs8176749", "rs8176720", "rs4751995", "rs10885997", "rs2303973", "rs75621879", "rs1784302", "rs540261", "rs1177562", "rs1047771", "rs2282133", "rs2296984", "rs3134586", "rs3134587", "rs3743248", "rs20543", "rs28584228", "rs7184958", "rs1869348", "rs2056822", "rs1744769", "rs6000174", "rs2227167", "rs12172195", "rs2142662", "rs1140555", "rs14136"};
         // ncbi.dbSNPIDESearch(snps);

          String httpsURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=((schizophrenia[tiab])+AND+(NM_004672[tiab]+OR+MAP3K6[tiab])+AND+(gene[tiab]+OR+genes[tiab]+OR+mRNA[tiab]+OR+protein[tiab]+OR+proteins[tiab]+OR+transcription[tiab]+OR+transcript[tiab]+OR+transcripts[tiab]+OR+expressed[tiab]+OR+expression[tiab]+OR+expressions[tiab]+OR+locus[tiab]+OR+loci[tiab]+OR+SNP[tiab]))&datetype=edat&retmax=100";
    URL myurl = new URL(httpsURL);
    HttpsURLConnection con = (HttpsURLConnection)myurl.openConnection();
    InputStream ins = con.getInputStream();
    InputStreamReader isr = new InputStreamReader(ins);
    BufferedReader in = new BufferedReader(isr);

    String inputLine;

    while ((inputLine = in.readLine()) != null)
    {
      System.out.println(inputLine);
    }

    in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
