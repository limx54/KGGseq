/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.List;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.hssf.util.HSSFColor;
import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.cobi.util.file.LocalFileFunc;

/**
 *
 * @author MX Li
 */
public class LocalExcelFile {

    /**
     * Creates a new instance of LocalExcelFile
     */
    public LocalExcelFile() {
    }

    public static void main(String[] args) {
        try {
            String inFile = "E:/home/mxli/MyJava/kggseq1/espEANoHM.flt.txt";
            String outFile = "E:/home/mxli/MyJava/kggseq1/espEANoHM.flt1.txt";
            BufferedReader br = new BufferedReader(new FileReader(inFile));
            BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
            String line = null;
            String[] cells = null;

            int[] orgIndices = new int[]{0, 1, 2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 33, 34, 35, 36, 37, 38, 39, 40};
            int selectedColNum = orgIndices.length;
            int i, pos;
            String delimiter = "\t";

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.trim().length() == 0) {
                    continue;
                }
                cells = line.split(delimiter, -1);
                if (cells[39].equals(".")) {
                    continue;
                }
                bw.write(cells[0]);
                for (i = 1; i < selectedColNum; i++) {
                    bw.write("\t");
                    bw.write(cells[orgIndices[i]]);
                }
                bw.write("\n");
            }
            bw.close();
            br.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public static boolean writeArray2ExcelSheet1(HSSFSheet sheet1, HSSFWorkbook wb, List<String[]> arry, boolean hasHead) throws
            Exception {
        int rowNum = arry.size();
        if (rowNum == 0) {
            System.err.println("No input data!");
            return false;
        }

        String[] titleNames = null;
        if (hasHead) {
            titleNames = (String[]) arry.get(0);
        }
        int columnNum = ((String[]) arry.get(0)).length;

        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth(i, (int) ((30 * 8) / ((double) 1 / 20)));
        }

        HSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        HSSFCellStyle headStyle = wb.createCellStyle();
        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(HSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        HSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        HSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.LIGHT_CORNFLOWER_BLUE.index);
        markedBodyStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        int rowIndex = 0;
        //create titile row
        HSSFRow row = sheet1.createRow(rowIndex);

        HSSFCell cell = null;
        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            rowIndex++;
        }

        for (int i = rowIndex; i < rowNum; i++) {
            row = sheet1.createRow(i);
            String[] line = (String[]) arry.get(i);
            columnNum = line.length;
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (line[0] != null) {
                    cell.setCellStyle(markedBodyStyle);
                } else {
                    cell.setCellStyle(bodyStyle);
                }
                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?hssf.?usermodel.?HSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue("");
                }

            }
        }

        return true;
    }

    public static boolean writeArray2ExcelFile1(String fileName, List<String[]> arry, boolean hasHead, int indexKey) throws
            Exception {
        int rowNum = arry.size();
        if (rowNum == 0) {
            System.err.println("No input data!");
            return false;
        }

        HSSFWorkbook wb = new HSSFWorkbook();
        HSSFSheet sheet1 = wb.createSheet("Data");
        String[] titleNames = null;
        if (hasHead) {
            titleNames = (String[]) arry.get(0);
        }
        int columnNum = ((String[]) arry.get(0)).length;

        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth(i, (int) ((30 * 6) / ((double) 1 / 20)));
        }

        HSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom font to the text in the comment
        HSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(HSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        HSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        HSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        int rowIndex = 0;
        //create titile row
        HSSFRow row = sheet1.createRow(rowIndex);

        String lastKey = null;
        int switcher = -1;
        HSSFCell cell = null;
        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            rowIndex++;
        }

        for (int i = rowIndex; i < rowNum; i++) {
            row = sheet1.createRow((i));
            String[] line = (String[]) arry.get(i);
            columnNum = line.length;
            if (indexKey >= 0) {
                if (lastKey == null && line[indexKey] != null) {
                    lastKey = line[indexKey];
                    switcher *= -1;
                } else if (lastKey != null && line[indexKey] == null) {
                    lastKey = line[indexKey];
                    switcher *= -1;
                } else if (lastKey == null && line[indexKey] == null) {
                } else {
                    if (!lastKey.equals(line[indexKey])) {
                        switcher *= -1;
                        lastKey = line[indexKey];
                    }
                }
            } else {
                switcher = 1;
            }
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (switcher > 0) {
                    cell.setCellStyle(bodyStyle);
                } else {
                    cell.setCellStyle(markedBodyStyle);
                }

                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?hssf.?usermodel.?HSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }

            }
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(fileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean writeMultArray2ExcelFile1(String fileName, List<List<String[]>> arrys, List<String> sheetLabels, boolean hasHead, int indexKey) throws
            Exception {

        if (arrys.isEmpty()) {
            System.err.println("No input data!");
            return false;
        }

        HSSFWorkbook wb = new HSSFWorkbook();
        HSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom font to the text in the comment
        HSSFFont font = wb.createFont();

        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(HSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        HSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);

        HSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(HSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle.setFillPattern(HSSFCellStyle.SOLID_FOREGROUND);

        String lastKey = null;
        int switcher = -1;
        HSSFCell cell = null;
        String[] titleNames = null;
        int d = 0;
        for (List<String[]> arry : arrys) {
            HSSFSheet sheet1 = wb.createSheet(sheetLabels.get(d));
            if (hasHead) {
                titleNames = (String[]) arry.get(0);
            }
            int columnNum = ((String[]) arry.get(0)).length;

            for (int i = 0; i < columnNum; i++) {
                sheet1.setColumnWidth(i, (int) ((30 * 6) / ((double) 1 / 20)));
            }

            int rowIndex = 0;
            //create titile row
            HSSFRow row = sheet1.createRow(rowIndex);

            if (titleNames != null) {
                for (int i = 0; i < columnNum; i++) {
                    cell = row.createCell(i);
                    cell.setCellValue(titleNames[i]);
                    cell.setCellStyle(headStyle);
                }
                rowIndex++;
            }
            int rowNum = arry.size();

            for (int i = rowIndex; i < rowNum; i++) {
                row = sheet1.createRow((i));
                String[] line = (String[]) arry.get(i);
                columnNum = line.length;
                if (indexKey >= 0) {
                    if (lastKey == null && line[indexKey] != null) {
                        lastKey = line[indexKey];
                        switcher *= -1;
                    } else if (lastKey != null && line[indexKey] == null) {
                        lastKey = line[indexKey];
                        switcher *= -1;
                    } else if (lastKey == null && line[indexKey] == null) {
                    } else {
                        if (!lastKey.equals(line[indexKey])) {
                            switcher *= -1;
                            lastKey = line[indexKey];
                        }
                    }
                } else {
                    switcher = 1;
                }
                for (int j = 0; j < columnNum; j++) {
                    cell = row.createCell(j);
                    if (switcher > 0) {
                        cell.setCellStyle(bodyStyle);
                    } else {
                        cell.setCellStyle(markedBodyStyle);
                    }

                    if (line[j] != null) {
                        if (Util.isNumeric(line[j])) {
                            //org.?apache.?poi.?hssf.?usermodel.?HSSFCell.CELL_TYPE_NUMERIC
                            cell.setCellType(0);
                            cell.setCellValue(Double.parseDouble(line[j]));
                        } else {
                            cell.setCellValue(line[j]);
                        }
                    } else {
                        cell.setCellValue(".");
                    }

                }
            }

            d++;
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(fileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean writeArray2ExcelFile1(String FileName, String[] titleNames, List<String[]> arry) throws
            Exception {
        HSSFWorkbook wb = new HSSFWorkbook();
        HSSFSheet sheet1 = wb.createSheet("Data");
        int columnNum = titleNames.length;
        int rowNum = arry.size();

        //for (int i = 0; i < columnNum; i++) {
        //  sheet1.setColumnWidth( i, (short) ((30 * 8) / ((double) 1 / 20)));
        //}
        HSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom headFont to the text in the comment
        HSSFFont headFont = wb.createFont();
        headFont.setFontName("Courier New");
        headFont.setFontHeightInPoints((short) 10);
        headFont.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
        headFont.setColor(HSSFColor.BLACK.index);

        headStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
        headStyle.setFillBackgroundColor(HSSFColor.GREY_25_PERCENT.index);
        headStyle.setFillForegroundColor(HSSFColor.BLACK.index);
        headStyle.setLocked(true);
        headStyle.setFont(headFont);
        headStyle.setBorderTop((short) 2);
        headStyle.setBorderBottom((short) 1);

        HSSFCellStyle contentStyle = wb.createCellStyle();
        //apply custom headFont to the text in the comment
        HSSFFont contentFont = wb.createFont();
        contentFont.setFontName("Courier New");
        contentFont.setFontHeightInPoints((short) 9);
        //headFont.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
        contentFont.setColor(HSSFColor.BLACK.index);

        contentStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
        contentStyle.setFillForegroundColor(HSSFColor.BLACK.index);
        contentStyle.setFont(contentFont);

        //create titile row
        HSSFRow row = sheet1.createRow(0);
        int heandLine = 0;
        HSSFCell cell = null;

        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            heandLine++;
        }

        for (int i = 0; i < rowNum; i++) {
            row = sheet1.createRow((i + heandLine));
            String[] line = (String[]) arry.get(i);
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?hssf.?usermodel.?HSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }

            }
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(FileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean writeArray2XLSXSheet(XSSFSheet sheet1, XSSFWorkbook wb, List<String[]> arry, boolean hasHead) throws
            Exception {
        int rowNum = arry.size();
        if (rowNum == 0) {
            System.err.println("No input data!");
            return false;
        }

        String[] titleNames = null;
        if (hasHead) {
            titleNames = (String[]) arry.get(0);
        }
        int columnNum = ((String[]) arry.get(0)).length;

        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth((short) i, (short) ((30 * 8) / ((double) 1 / 20)));
        }

        XSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        XSSFCellStyle headStyle = wb.createCellStyle();
        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(XSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        XSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        XSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.LIGHT_CORNFLOWER_BLUE.index);
        markedBodyStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        int rowIndex = 0;
        //create titile row
        XSSFRow row = sheet1.createRow(rowIndex);

        XSSFCell cell = null;
        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell((short) i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            rowIndex++;
        }

        for (int i = rowIndex; i < rowNum; i++) {
            row = sheet1.createRow((i));
            String[] line = (String[]) arry.get(i);
            columnNum = line.length;
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (line[0] != null) {
                    cell.setCellStyle(markedBodyStyle);
                } else {
                    cell.setCellStyle(bodyStyle);
                }
                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue("");
                }

            }
        }
        return true;
    }

    public static boolean writeArray2XLSXFile(String fileName, List<String[]> arry, boolean hasHead, int indexKey1, int indexKey2) throws
            Exception {
        int rowNum = arry.size();
        if (rowNum == 0) {
            System.err.println("No input data!");
            return false;
        }

        XSSFWorkbook wb = new XSSFWorkbook();
        XSSFSheet sheet1 = wb.createSheet("Data");
        String[] titleNames = null;
        if (hasHead) {
            titleNames = (String[]) arry.get(0);
        }
        int columnNum = ((String[]) arry.get(0)).length;

        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth(i, (short) ((30 * 6) / ((double) 1 / 20)));
        }

        XSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom font to the text in the comment
        XSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(XSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        XSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        XSSFCellStyle markedBodyStyle1 = wb.createCellStyle();
        markedBodyStyle1.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle1.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle1.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle1.setFillForegroundColor(HSSFColor.LIGHT_YELLOW.index);
        markedBodyStyle1.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        XSSFCellStyle markedBodyStyle2 = wb.createCellStyle();
        markedBodyStyle2.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle2.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle2.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle2.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle2.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        int rowIndex = 0;
        //create titile row
        XSSFRow row = sheet1.createRow(rowIndex);

        String lastKey = null;
        int switcher = -1;
        XSSFCell cell = null;
        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            rowIndex++;
        }

        for (int i = rowIndex; i < rowNum; i++) {
            row = sheet1.createRow((i));
            String[] line = arry.get(i);

            columnNum = line.length;
            if (indexKey2 >= 0) {
                if (lastKey == null && line[indexKey2] != null) {
                    lastKey = line[indexKey2];
                    switcher *= -1;
                } else if (lastKey != null && line[indexKey2] == null) {
                    lastKey = line[indexKey2];
                    switcher *= -1;
                } else if (lastKey == null && line[indexKey2] == null) {
                } else {
                    if (!lastKey.equals(line[indexKey2])) {
                        switcher *= -1;
                        lastKey = line[indexKey2];
                    }
                }
            } else {
                switcher = 1;
            }
            // System.out.println(cells1[0]);
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (indexKey1 >= 0 && line[indexKey1] != null) {
                    cell.setCellStyle(markedBodyStyle1);
                } else {
                    if (switcher > 0) {
                        cell.setCellStyle(bodyStyle);
                    } else {
                        cell.setCellStyle(markedBodyStyle2);
                    }
                }

                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }

            }
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(fileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean convertTextFile2XLSXFile(String inFileName, String outFileName, boolean hasHead, int indexKey) throws
            Exception {
        BufferedReader br = LocalFileFunc.getBufferedReader(inFileName);
        String line = br.readLine();
        if (line == null) {
            return false;
        }
        String[] cells1 = Util.tokenize(line, '\t');
        XSSFWorkbook wb = new XSSFWorkbook();
        XSSFSheet sheet1 = wb.createSheet("Data");

        int columnNum = cells1.length;
        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth(i, (short) ((30 * 6) / ((double) 1 / 20)));
        }

        XSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom font to the text in the comment
        XSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(XSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        XSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        XSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        int rowIndex = 0;
        //create titile row
        XSSFRow row = sheet1.createRow(rowIndex);

        String lastKey = null;
        int switcher = -1;
        XSSFCell cell = null;
        if (hasHead) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(cells1[i]);
                cell.setCellStyle(headStyle);
            }
        } else {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(cells1[i]);
                cell.setCellStyle(bodyStyle);
            }
        }
        rowIndex++;
        while ((line = br.readLine()) != null) {
            cells1 = Util.tokenize(line, '\t');
            row = sheet1.createRow((rowIndex));
            columnNum = cells1.length;
            if (indexKey >= 0) {
                if (lastKey == null && cells1[indexKey] != null) {
                    lastKey = cells1[indexKey];
                    switcher *= -1;
                } else if (lastKey != null && cells1[indexKey] == null) {
                    lastKey = cells1[indexKey];
                    switcher *= -1;
                } else if (lastKey == null && cells1[indexKey] == null) {
                } else {
                    if (!lastKey.equals(cells1[indexKey])) {
                        switcher *= -1;
                        lastKey = cells1[indexKey];
                    }
                }
            } else {
                switcher = 1;
            }
            // System.out.println(cells1[0]);
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (switcher > 0) {
                    cell.setCellStyle(bodyStyle);
                } else {
                    cell.setCellStyle(markedBodyStyle);
                }

                if (cells1[j] != null) {
                    if (Util.isNumeric(cells1[j])) {
                        //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(cells1[j]));
                    } else {
                        cell.setCellValue(cells1[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }
            }
            rowIndex++;
        }
        br.close();
        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(outFileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean appendArray2XLSXFile(String fileName, List<String[]> arry, boolean hasHead, int indexKey) throws
            Exception {
        int rowNum = arry.size();
        if (rowNum == 0) {
            System.err.println("No input data!");
            return false;
        }
        XSSFWorkbook wb;
        XSSFSheet sheet1;
        File file = new File(fileName);
        if (!file.exists()) {
            file.getParentFile().mkdirs();
            wb = new XSSFWorkbook();
            sheet1 = wb.createSheet("Data");
        } else {
            if (hasHead) {
                //I hope this is a new file this time
                file.delete();
                wb = new XSSFWorkbook();
                sheet1 = wb.createSheet("Data");
            } else {
                FileInputStream inFile = new FileInputStream(fileName);
                wb = new XSSFWorkbook(inFile);
                sheet1 = wb.getSheetAt(0);
                inFile.close();
            }

        }

        String[] titleNames = null;
        int columnNum = ((String[]) arry.get(0)).length;
        for (int i = 0; i < columnNum; i++) {
            sheet1.setColumnWidth(i, (short) ((30 * 6) / ((double) 1 / 20)));
        }

        //apply custom font to the text in the comment
        XSSFFont font = wb.createFont();
        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        XSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        XSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        int overallRowIndex = 0;
        if (!hasHead) {
            overallRowIndex = sheet1.getLastRowNum() + 1;
            rowNum = rowNum + overallRowIndex;
        }

        int currRowIndex = 0;
        //create titile row
        XSSFRow row = sheet1.createRow(overallRowIndex);

        String lastKey = null;
        int switcher = -1;
        XSSFCell cell = null;
        if (hasHead && titleNames != null) {
            titleNames = (String[]) arry.get(0);
            XSSFCellStyle headStyle = wb.createCellStyle();
            headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
            headStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

            headStyle.setBorderBottom(XSSFCellStyle.BORDER_THIN);
            headStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
            headStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
            headStyle.setBorderTop(XSSFCellStyle.BORDER_THIN);
            headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

            headStyle.setLocked(true);
            headStyle.setFont(font);
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            overallRowIndex++;
            currRowIndex++;
        }

        for (int i = overallRowIndex; i < rowNum; i++) {
            row = sheet1.createRow((i));
            String[] line = arry.get(currRowIndex);

            columnNum = line.length;
            if (indexKey >= 0) {
                if (lastKey == null && line[indexKey] != null) {
                    lastKey = line[indexKey];
                    switcher *= -1;
                } else if (lastKey != null && line[indexKey] == null) {
                    lastKey = line[indexKey];
                    switcher *= -1;
                } else if (lastKey == null && line[indexKey] == null) {
                } else {
                    if (!lastKey.equals(line[indexKey])) {
                        switcher *= -1;
                        lastKey = line[indexKey];
                    }
                }
            } else {
                switcher = 1;
            }
            // System.out.println(cells1[0]);
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (switcher > 0) {
                    cell.setCellStyle(bodyStyle);
                } else {
                    cell.setCellStyle(markedBodyStyle);
                }

                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }

            }
            currRowIndex++;
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(fileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean writeMultArray2XLSXFile(String fileName, List<List<String[]>> arrys, List<String> sheetLabels, boolean hasHead, int indexKey) throws
            Exception {

        if (arrys.isEmpty()) {
            System.err.println("No input data!");
            return false;
        }

        XSSFWorkbook wb = new XSSFWorkbook();
        XSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom font to the text in the comment
        XSSFFont font = wb.createFont();

        font.setFontName("Courier New");
        font.setFontHeightInPoints((short) 10);
        font.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        font.setColor(HSSFColor.RED.index);

        headStyle.setFillForegroundColor(HSSFColor.LIGHT_GREEN.index);
        headStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        headStyle.setBorderBottom(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        headStyle.setBorderTop(XSSFCellStyle.BORDER_THIN);
        headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        headStyle.setLocked(true);
        headStyle.setFont(font);

        XSSFCellStyle bodyStyle = wb.createCellStyle();
        bodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        bodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);

        XSSFCellStyle markedBodyStyle = wb.createCellStyle();
        markedBodyStyle.setBorderLeft(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setBorderRight(XSSFCellStyle.BORDER_THIN);
        markedBodyStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        markedBodyStyle.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
        markedBodyStyle.setFillPattern(XSSFCellStyle.SOLID_FOREGROUND);

        String lastKey = null;
        int switcher = -1;
        XSSFCell cell = null;
        String[] titleNames = null;
        int d = 0;
        for (List<String[]> arry : arrys) {
            XSSFSheet sheet1 = wb.createSheet(sheetLabels.get(d));
            if (hasHead) {
                titleNames = (String[]) arry.get(0);
            }
            int columnNum = ((String[]) arry.get(0)).length;

            for (int i = 0; i < columnNum; i++) {
                sheet1.setColumnWidth(i, (short) ((30 * 6) / ((double) 1 / 20)));
            }

            int rowIndex = 0;
            //create titile row
            XSSFRow row = sheet1.createRow(rowIndex);

            if (titleNames != null) {
                for (int i = 0; i < columnNum; i++) {
                    cell = row.createCell(i);
                    cell.setCellValue(titleNames[i]);
                    cell.setCellStyle(headStyle);
                }
                rowIndex++;
            }
            int rowNum = arry.size();

            for (int i = rowIndex; i < rowNum; i++) {
                row = sheet1.createRow((i));
                String[] line = (String[]) arry.get(i);
                columnNum = line.length;
                if (indexKey >= 0) {
                    if (lastKey == null && line[indexKey] != null) {
                        lastKey = line[indexKey];
                        switcher *= -1;
                    } else if (lastKey != null && line[indexKey] == null) {
                        lastKey = line[indexKey];
                        switcher *= -1;
                    } else if (lastKey == null && line[indexKey] == null) {
                    } else {
                        if (!lastKey.equals(line[indexKey])) {
                            switcher *= -1;
                            lastKey = line[indexKey];
                        }
                    }
                } else {
                    switcher = 1;
                }
                for (int j = 0; j < columnNum; j++) {
                    cell = row.createCell(j);
                    if (switcher > 0) {
                        cell.setCellStyle(bodyStyle);
                    } else {
                        cell.setCellStyle(markedBodyStyle);
                    }

                    if (line[j] != null) {
                        if (Util.isNumeric(line[j])) {
                            //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                            cell.setCellType(0);
                            cell.setCellValue(Double.parseDouble(line[j]));
                        } else {
                            cell.setCellValue(line[j]);
                        }
                    } else {
                        cell.setCellValue(".");
                    }

                }
            }

            d++;
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(fileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

    public static boolean writeArray2XLSXFile(String FileName, String[] titleNames, List<String[]> arry) throws
            Exception {
        XSSFWorkbook wb = new XSSFWorkbook();
        XSSFSheet sheet1 = wb.createSheet("Data");
        int columnNum = titleNames.length;
        int rowNum = arry.size();

        //for (int i = 0; i < columnNum; i++) {
        //  sheet1.setColumnWidth( i, (short) ((30 * 8) / ((double) 1 / 20)));
        //}
        XSSFCellStyle headStyle = wb.createCellStyle();
        //apply custom headFont to the text in the comment
        XSSFFont headFont = wb.createFont();
        headFont.setFontName("Courier New");
        headFont.setFontHeightInPoints((short) 10);
        headFont.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        headFont.setColor(HSSFColor.BLACK.index);

        headStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        headStyle.setFillBackgroundColor(HSSFColor.GREY_25_PERCENT.index);
        headStyle.setFillForegroundColor(HSSFColor.BLACK.index);
        headStyle.setLocked(true);
        headStyle.setFont(headFont);
        headStyle.setBorderTop((short) 2);
        headStyle.setBorderBottom((short) 1);

        XSSFCellStyle contentStyle = wb.createCellStyle();
        //apply custom headFont to the text in the comment
        XSSFFont contentFont = wb.createFont();
        contentFont.setFontName("Courier New");
        contentFont.setFontHeightInPoints((short) 9);
        //headFont.setBoldweight(XSSFFont.BOLDWEIGHT_BOLD);
        contentFont.setColor(HSSFColor.BLACK.index);

        contentStyle.setAlignment(XSSFCellStyle.ALIGN_CENTER);
        contentStyle.setFillForegroundColor(HSSFColor.BLACK.index);
        contentStyle.setFont(contentFont);

        //create titile row
        XSSFRow row = sheet1.createRow(0);
        int heandLine = 0;
        XSSFCell cell = null;

        if (titleNames != null) {
            for (int i = 0; i < columnNum; i++) {
                cell = row.createCell(i);
                cell.setCellValue(titleNames[i]);
                cell.setCellStyle(headStyle);
            }
            heandLine++;
        }

        for (int i = 0; i < rowNum; i++) {
            row = sheet1.createRow((i + heandLine));
            String[] line = (String[]) arry.get(i);
            for (int j = 0; j < columnNum; j++) {
                cell = row.createCell(j);
                if (line[j] != null) {
                    if (Util.isNumeric(line[j])) {
                        //org.?apache.?poi.?XSSF.?usermodel.?XSSFCell.CELL_TYPE_NUMERIC
                        cell.setCellType(0);
                        cell.setCellValue(Double.parseDouble(line[j]));
                    } else {
                        cell.setCellValue(line[j]);
                    }
                } else {
                    cell.setCellValue(".");
                }

            }
        }

        // Write the output to a inFile
        FileOutputStream fileOut = new FileOutputStream(FileName);
        wb.write(fileOut);
        fileOut.close();

        return true;
    }

}
