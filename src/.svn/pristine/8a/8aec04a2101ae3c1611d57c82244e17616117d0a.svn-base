// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
package org.cobi.util.plot;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.stat.Descriptive;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;
import umontreal.iro.lecuyer.probdist.BetaDist;

/**
 *
 * @author mxli
 */
public class PValuePainter extends BasicPainter {

    public PValuePainter(int width, int height) {
        super(width, height);
    }

    /**
     *
     * @param valueList
     * @param title
     * @param outputPath
     * @param pValueTolerationLevle
     * @return
     * @throws Exception
     */
    public int drawQQPlot(final DoubleArrayList valueList, String title, String outputPath, double pValueTolerationLevle) throws Exception {
        if (valueList == null || valueList.size() == 0) {
            System.err.println("Null p-value list");
            return -1;
        }
        if (title != null && title.length() > 0) {
            calculateDataPlottingArea(true);
        } else {
            calculateDataPlottingArea(false);
        }
        DoubleArrayList tmpValueList = valueList.copy();
        BufferedImage image = new BufferedImage(this.canvasWidth, this.canvasHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintCanvas(g2d);


        int dataSize = tmpValueList.size();
        tmpValueList.quickSort();
        //max and min value in vertical
        double vMax = -Math.log10(tmpValueList.getQuick(0));
        double vMin = -Math.log10(tmpValueList.getQuick(dataSize - 1));
        double hMin = 0;
        double hMax = -Math.log10((1.0 / (dataSize + 1)));

        //sometimes the p values are too large
        if (vMax > -Math.log10(pValueTolerationLevle)) {
            vMax = -Math.log10(pValueTolerationLevle);
        }

        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
                dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
        int maxYScalLen = drawAxesScale(g2d, transformer);
        drawAxes(g2d, title, "Expected [-log10(P)]", "Observed [-log10(P)]", maxYScalLen);
        //draw the standard line
        Point point1 = transformer.data2ScreenPoint(0, 0);
        Point point2 = transformer.data2ScreenPoint(Math.min(hMax, vMax), Math.min(hMax, vMax));
        Line2D.Double zz = new Line2D.Double(point1, point2);
        Stroke oldStroke = g2d.getStroke();
        Color oldColor = g2d.getColor();
        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2d.draw(zz);
        g2d.setStroke(oldStroke);
        g2d.setColor(oldColor);

        int rectangleSize = 2;
        double maxDiff = 0.05;
        int maxDiffIndex = -1;
        double tmpDouble2;
        double tmpDouble1;
        double triangleLen = 4;
        int dataSizeMore = dataSize + 1;
        for (int i = 0; i < dataSize; i++) {
            tmpDouble1 = -Math.log10((double) (i + 1) / (dataSizeMore));
            tmpDouble2 = -Math.log10(tmpValueList.getQuick(i));
            if (tmpDouble2 <= vMax) {
                transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                g2d.drawRect(point1.x - rectangleSize / 2, point1.y - rectangleSize / 2, rectangleSize, rectangleSize);
                //Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,  point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                // g2d.draw(f);
            } else {
                tmpDouble2 = vMax;
                transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                double x = point1.getX();
                double y = point1.getY();

                g2d.setColor(Color.red);
                plotTriangle(g2d, (int) x, (int) y, (int) triangleLen);
                g2d.setColor(oldColor);
            }


            /*
            if (maxDiffIndex < 0 && (Math.abs(expectedList.getQuick(i) - tmpValueList.getQuick(i)) > maxDiff)) {
            maxDiffIndex = i;
            }
             *
             */
        }

        /*
        float[] dashes = {3.f};
        g2d.setColor(Color.RED);
        g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
        
        //draw large difference
        point1 = transformer.data2ScreenPoint(expectedList.getQuick(maxDiffIndex), tmpValueList.getQuick(maxDiffIndex));
        point2 = transformer.data2ScreenPoint(hMax, vMax);
        double vRectangleSize = point1.getY() - point2.getY();
        double hRectangleSize = point2.getX() - point1.getX();
        //zz = new Line2D.Double(point1, point2);
        // g2d.draw(zz);
        Rectangle2D.Double f = new Rectangle2D.Double(point1.getX(),
        point1.getY() - vRectangleSize, hRectangleSize, vRectangleSize);
        g2d.draw(f);
        
        g2d.setStroke(oldStroke);
        g2d.setColor(oldColor);
         */

        g2d.dispose();
        outputPNGFile(image, outputPath);
        return maxDiffIndex;
    }

    /**
     *Do QQ plot by JFreeChart. But it is very very slow for large dataset.
     * @param valueList
     * @param title
     * @param outputPath
     * @param pValueTolerationLevle
     * @return
     * @throws Exception
     */
    /*
    public int drawQQPlotJFreeChart(final DoubleArrayList valueList, String title, String outputPath, double pValueTolerationLevle) throws Exception {
    if (valueList == null || valueList.size() == 0) {
    System.err.println("Null p-value list");
    return -1;
    }
    XYSeriesCollection dataset = new XYSeriesCollection();
    XYSeries series = new XYSeries("");
    
    DoubleArrayList tmpValueList = valueList.copy();
    int dataSize = tmpValueList.size();
    tmpValueList.quickSort();
    //max and min value in vertical
    double vMax = -Math.log10(tmpValueList.getQuick(0));
    double vMin = -Math.log10(tmpValueList.getQuick(dataSize - 1));
    double hMin = 0;
    double hMax = -Math.log10((1.0 / (dataSize + 1)));
    //sometimes the p values are too large
    if (vMax > -Math.log10(pValueTolerationLevle)) {
    vMax = -Math.log10(pValueTolerationLevle);
    }
    double tmpDouble2;
    double tmpDouble1;
    int dataSizeMore = dataSize + 1;
    int dataSizeLess = dataSize - 1;
    for (int i = 0; i < dataSizeLess; i++) {
    tmpDouble1 = -Math.log10((double) (i + 1) / (dataSizeMore));
    tmpDouble2 = -Math.log10(tmpValueList.getQuick(i));
    if (tmpDouble2 <= vMax) {
    series.add(tmpDouble1, tmpDouble2, false);
    } else {
    tmpDouble2 = vMax;
    series.add(tmpDouble1, tmpDouble2, false);
    }
    }
    tmpDouble1 = -Math.log10((double) (dataSizeLess + 1) / (dataSizeMore));
    tmpDouble2 = -Math.log10(tmpValueList.getQuick(dataSizeLess));
    if (tmpDouble2 <= vMax) {
    series.add(tmpDouble1, tmpDouble2);
    } else {
    tmpDouble2 = vMax;
    series.add(tmpDouble1, tmpDouble2);
    }
    dataset.addSeries(series);
    JFreeChart chart = ChartFactory.createXYLineChart(title, "Expected (-log(10)P)", "Observed (-log(10)P)", dataset, PlotOrientation.VERTICAL, false, false, false);
    
    BufferedImage image = chart.createBufferedImage(this.figureWidth, this.figureHeight);
    
    
    Graphics2D g2d = image.createGraphics();
    CoordinateTransformer transformer = new CoordinateTransformer();
    transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
    dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
    
    //draw the standard line
    Point2D point1 = transformer.data2ScreenPoint(0, 0);
    Point2D point2 = transformer.data2ScreenPoint(Math.max(hMax, vMax), Math.max(hMax, vMax));
    Line2D.Double zz = new Line2D.Double(point1, point2);
    Stroke oldStroke = g2d.getStroke();
    Color oldColor = g2d.getColor();
    g2d.setColor(Color.GREEN);
    g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
    g2d.draw(zz);
    g2d.setStroke(oldStroke);
    g2d.setColor(oldColor);
    
    g2d.dispose();
    //outputJPEGFile(image, outputPath);
    outputPNGFile(image, outputPath);
    return 0;
    }
     */
    /**
     *
     * @param valueLists
     * @param legends
     * @param title
     * @param outputPath
     * @param pValueTolerationLevle
     * @return
     * @throws Exception
     */
    public void drawMultipleQQPlot(final List<DoubleArrayList> valueLists, List<String> legends,
            String title, String outputPath, double pValueTolerationLevle) throws Exception {
        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return;
        }

        int listNum = valueLists.size();
        for (int i = listNum - 1; i >= 0; i--) {
            if (valueLists.get(i) == null || valueLists.get(i).size() == 0) {
                System.err.println("Null p-value list");
                valueLists.remove(i);
            }
        }

        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return;
        }

        int dataSize = valueLists.get(0).size();
        if (title != null && title.length() > 0) {
            calculateDataPlottingArea(true);
        } else {
            calculateDataPlottingArea(false);
        }
        //max and min value in vertical
        double vMax = (valueLists.get(0).getQuick(0));
        double vMin = vMax;
        double hMin = 0;
        double hMax = (1.0 / dataSize);

        int maxDataPointSize = dataSize;
        for (int i = 0; i < listNum; i++) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            tmpValueList.quickSort();
            dataSize = tmpValueList.size();
            for (int j = 0; j < dataSize; j++) {
                if (vMax > tmpValueList.getQuick(j)) {
                    vMax = tmpValueList.getQuick(j);
                } else if (vMin < tmpValueList.getQuick(j)) {
                    //vMin = tmpValueList.getQuick(j);
                }
            }
            if (hMax > (1.0 / (dataSize + 1))) {
                hMax = (1.0 / (dataSize + 1));
            }
            if (maxDataPointSize < dataSize) {
                maxDataPointSize = dataSize;
            }
        }

        vMin=1;
        //convert them into the data to be used
        vMin = -Math.log10(vMin);
        vMax = -Math.log10(vMax);
        hMax = -Math.log10(hMax);

        //sometimes the p values are too large
        if (vMax > -Math.log10(pValueTolerationLevle)) {
            vMax = -Math.log10(pValueTolerationLevle);
        }

        BufferedImage image = new BufferedImage(canvasWidth, canvasHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintCanvas(g2d);

        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
                dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
        int maxYScalLen = drawAxesScale(g2d, transformer);
        drawAxes(g2d, title, "Expected [-log10(P)]", "Observed [-log10(P)]", maxYScalLen);
        //draw the standard line
        Point point1 = transformer.data2ScreenPoint(0, 0);
        Point point2 = transformer.data2ScreenPoint(Math.min(hMax, vMax), Math.min(hMax, vMax));

        Line2D.Double zz = new Line2D.Double(point1, point2);
        Stroke oldStroke = g2d.getStroke();
        Color oldColor = g2d.getColor();
        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2d.draw(zz);

        int rectangleSize = 2;
        int harlfRectangleSize = rectangleSize / 2;
        int triangleLen = 4;


        double maxDiff = 0.05;
        Color[] colors = {Color.BLUE, Color.RED, Color.ORANGE, Color.MAGENTA, Color.BLACK, Color.PINK};

        //Color[] colors = {new Color(00, 00, 0xff), new Color(0xff, 00, 00), new Color(00, 0x88, 0xff),    new Color(0xff, 0x88, 00), new Color(0x88, 00, 0xff), new Color(00, 00, 00)};
        int strHei = g2d.getFontMetrics().getHeight();
        double vRectangleSize = 6;
        double hRectangleSize = 6;
        double tmpDouble1;
        double tmpDouble2;
        g2d.setStroke(oldStroke);

        Font oldFront = g2d.getFont();
        g2d.setFont(numFont);
        int maxLegentWidth = g2d.getFontMetrics().stringWidth(legends.get(listNum - 1));
        for (int i = listNum - 2; i >= 0; i--) {
            int strLen = g2d.getFontMetrics().stringWidth(legends.get(i));
            if (maxLegentWidth < strLen) {
                maxLegentWidth = strLen;
            }
        }

        double x, fmin, fmax;
        BetaDist betaDist = null;

        Point minP1 = new Point(), minP2 = new Point();
        Point maxP1 = new Point(), maxP2 = new Point();

        betaDist = new BetaDist(1, maxDataPointSize);
        x = betaDist.inverseF(0.5);

        x = -Math.log10(x);
        fmin = betaDist.inverseF(.025);
        fmin = -Math.log10(fmin);
        transformer.data2ScreenPoint(minP1, x, fmin);

        fmax = betaDist.inverseF(.975);
        fmax = -Math.log10(fmax);
        transformer.data2ScreenPoint(maxP1, x, fmax);

        g2d.setColor(Color.BLACK);
        float[] dashes = {3.f};
        g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
        // maxDataPointSize
        for (int k = 2; k <= maxDataPointSize;) {
            betaDist = new BetaDist(k, maxDataPointSize + 1 - k);
            x = betaDist.inverseF(0.5);
            x = -Math.log10(x);
            fmin = betaDist.inverseF(.025);
            fmin = -Math.log10(fmin);

            transformer.data2ScreenPoint(minP2, x, fmin);
            g2d.drawLine(minP1.x, minP1.y, minP2.x, minP2.y);
            minP1.x = minP2.x;
            minP1.y = minP2.y;
            fmax = betaDist.inverseF(.975);
            fmax = -Math.log10(fmax);

            transformer.data2ScreenPoint(maxP2, x, fmax);
            g2d.drawLine(maxP1.x, maxP1.y, maxP2.x, maxP2.y);

            maxP1.x = maxP2.x;
            maxP1.y = maxP2.y;
            if (k <= 2000) {
                k++;
            } else {
                k += k;
            }
        }
        g2d.setStroke(oldStroke);
        for (int i = listNum - 1; i >= 0; i--) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            dataSize = tmpValueList.size();
            g2d.setColor(colors[i]);
            int dataSizeMore = dataSize + 1;
            for (int j = 0; j < dataSize; j++) {
                tmpDouble1 = -Math.log10((double) (j + 1) / (dataSizeMore));
                tmpDouble2 = -Math.log10(tmpValueList.getQuick(j));
                if (tmpDouble2 <= vMax) {
                    transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                    g2d.drawRect(point1.x - harlfRectangleSize, point1.y - harlfRectangleSize, rectangleSize, rectangleSize);
                    //Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,  point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                    // g2d.draw(f);
                } else {
                    tmpDouble2 = vMax;
                    transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                    x = point1.getX();
                    double y = point1.getY();
                    plotTriangle(g2d, (int) x, (int) y, (int) triangleLen);
                }
            }
            //plot legends

            oldStroke = g2d.getStroke();

            //g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
            g2d.setStroke(new BasicStroke(1.5f));
            double xOffset = dataPlottingArea.getWidth() / 4;
            double yOffset = dataPlottingArea.getHeight() / 8;
            g2d.drawString(legends.get(i), (float) (dataPlottingArea.x + xOffset),
                    (float) (dataPlottingArea.getY() + yOffset) + strHei * i);
            Rectangle2D.Double f = new Rectangle2D.Double(dataPlottingArea.x + xOffset + maxLegentWidth + hRectangleSize,
                    dataPlottingArea.getY() + yOffset + strHei * i - vRectangleSize, hRectangleSize, vRectangleSize);

            g2d.fill(f);
            g2d.draw(f);
            g2d.setStroke(oldStroke);

        }
        g2d.setFont(oldFront);
        g2d.setColor(oldColor);

        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);

    }

    public static void main(String[] args) {
        try {

            DoubleArrayList dlist = new DoubleArrayList();
            Uniform twister = new Uniform(new MersenneTwister(new java.util.Date()));
            for (int i = 0; i < 1000000; ++i) {
                dlist.add(twister.nextDouble());
            }
            PValuePainter painter = new PValuePainter(400, 400);
            double pValueTolerationLevle = 1E-20;
            List<String> legs = new ArrayList<String>();
            legs.add("test");
            List<DoubleArrayList> pList = new ArrayList<DoubleArrayList>();
            pList.add(dlist);
            painter.drawMultipleQQPlot(pList, legs, null, "qq.png", pValueTolerationLevle);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void drawMultipleTruncatedQQPlot(final List<DoubleArrayList> valueLists, List<String> legends,
            String title, String outputPath, double pValueTolerationLevle, int maxPointNum) throws Exception {
        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return;
        }

        int listNum = valueLists.size();
        for (int i = listNum - 1; i >= 0; i--) {
            if (valueLists.get(i) == null || valueLists.get(i).size() == 0) {
                System.err.println("Null p-value list");
                valueLists.remove(i);
            }
        }

        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return;
        }

        int dataSize = valueLists.get(0).size();
        if (title != null && title.length() > 0) {
            calculateDataPlottingArea(true);
        } else {
            calculateDataPlottingArea(false);
        }
        //max and min value in vertical
        double vMax = (valueLists.get(0).getQuick(0));
        double vMin = vMax;
        double hMin = 0;
        double hMax = (1.0 / dataSize);
        double minMedian = Double.MAX_VALUE;
        int minMedianIndex = 0;
        double tmpDouble = 0;

        for (int i = 0; i < listNum; i++) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            tmpValueList.quickSort();
            tmpDouble = Descriptive.median(tmpValueList);
            dataSize = tmpValueList.size();
            if (minMedian > tmpDouble) {
                minMedian = tmpDouble;
                minMedianIndex = tmpValueList.binarySearch(minMedian);
                if (minMedianIndex < 0) {
                    minMedianIndex = -minMedianIndex - 1;
                }
            }

            for (int j = 0; j < dataSize; j++) {
                if (vMax > tmpValueList.getQuick(j)) {
                    vMax = tmpValueList.getQuick(j);
                } else if (vMin < tmpValueList.getQuick(j)) {
                    vMin = tmpValueList.getQuick(j);
                }
            }
            if (hMax > (1.0 / (dataSize + 1))) {
                hMax = (1.0 / (dataSize + 1));
            }
        }

        if (hMax > (1.0 / (maxPointNum + 1))) {
            hMax = (1.0 / (maxPointNum + 1));
        }

        //int sizeDiff = (int) (1 / minMedian) - minMedianIndex;
        int sizeDiff = 4000;
        //convert them into the data to be used
        //vMin = -Math.log10(vMin);
        vMin = 0;
        vMax = -Math.log10(vMax);
        hMax = -Math.log10(hMax);

        //sometimes the p values are too large
        if (vMax > -Math.log10(pValueTolerationLevle)) {
            vMax = -Math.log10(pValueTolerationLevle);
        }

        BufferedImage image = new BufferedImage(canvasWidth, canvasHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintCanvas(g2d);

        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
                dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
        int maxYScalLen = drawAxesScale(g2d, transformer);
        drawAxes(g2d, title, "Expected [-log10(P)]", "Observed [-log10(P)]", maxYScalLen);
        //draw the standard line
        Point point1 = transformer.data2ScreenPoint(0, 0);
        Point point2 = transformer.data2ScreenPoint(Math.min(hMax, vMax), Math.min(hMax, vMax));

        Line2D.Double zz = new Line2D.Double(point1, point2);
        Stroke oldStroke = g2d.getStroke();
        Color oldColor = g2d.getColor();
        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2d.draw(zz);

        int rectangleSize = 2;
        int harlfRectangleSize = rectangleSize / 2;
        int triangleLen = 4;


        double maxDiff = 0.05;
        Color[] colors = {Color.BLUE, Color.RED, Color.ORANGE, Color.MAGENTA, Color.BLACK, Color.PINK};

        //Color[] colors = {new Color(00, 00, 0xff), new Color(0xff, 00, 00), new Color(00, 0x88, 0xff),    new Color(0xff, 0x88, 00), new Color(0x88, 00, 0xff), new Color(00, 00, 00)};
        int strHei = g2d.getFontMetrics().getHeight();
        double vRectangleSize = 6;
        double hRectangleSize = 6;
        double tmpDouble1;
        double tmpDouble2;
        g2d.setStroke(oldStroke);

        Font oldFront = g2d.getFont();
        g2d.setFont(numFont);
        int maxLegentWidth = g2d.getFontMetrics().stringWidth(legends.get(listNum - 1));
        for (int i = listNum - 2; i >= 0; i--) {
            int strLen = g2d.getFontMetrics().stringWidth(legends.get(i));
            if (maxLegentWidth < strLen) {
                maxLegentWidth = strLen;
            }
        }

        for (int i = listNum - 1; i >= 0; i--) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            dataSize = tmpValueList.size();

            g2d.setColor(colors[i]);
            int dataSizeMore = dataSize + 1;
            for (int j = 0; j < dataSize; j++) {
                tmpDouble1 = -Math.log10((double) (j + 1) / (dataSizeMore + sizeDiff));
                tmpDouble2 = -Math.log10(tmpValueList.getQuick(j));
                if (tmpDouble2 <= vMax) {
                    transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                    g2d.drawRect(point1.x - harlfRectangleSize, point1.y - harlfRectangleSize, rectangleSize, rectangleSize);
                    //Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,  point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                    // g2d.draw(f);
                } else {
                    tmpDouble2 = vMax;
                    transformer.data2ScreenPoint(point1, tmpDouble1, tmpDouble2);
                    double x = point1.getX();
                    double y = point1.getY();
                    plotTriangle(g2d, (int) x, (int) y, (int) triangleLen);
                }
            }
            //plot legends
            float[] dashes = {3.f};
            oldStroke = g2d.getStroke();

            //g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
            g2d.setStroke(new BasicStroke(1.5f));
            double xOffset = dataPlottingArea.getWidth() / 4;
            double yOffset = dataPlottingArea.getHeight() / 8;
            g2d.drawString(legends.get(i), (float) (dataPlottingArea.x + xOffset),
                    (float) (dataPlottingArea.getY() + yOffset) + strHei * i);
            Rectangle2D.Double f = new Rectangle2D.Double(dataPlottingArea.x + xOffset + maxLegentWidth + hRectangleSize,
                    dataPlottingArea.getY() + yOffset + strHei * i - vRectangleSize, hRectangleSize, vRectangleSize);

            g2d.fill(f);
            g2d.draw(f);
            g2d.setStroke(oldStroke);

        }
        g2d.setFont(oldFront);
        g2d.setColor(oldColor);

        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);

    }

    /**
     *
     * @param xList
     * @param yList
     * @param title
     * @param outputPath
     * @throws Exception
     */
    public void drawScatterPlot(DoubleArrayList xList, DoubleArrayList yList, String title, String outputPath) throws Exception {
        BufferedImage image = new BufferedImage(this.canvasWidth, this.canvasHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintCanvas(g2d);


        g2d.setStroke(new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        int pointSize = xList.size();
        //max and min value in vertical
        double vMax = Descriptive.max(yList);
        double vMin = Descriptive.min(yList);
        double hMin = Descriptive.min(xList);
        double hMax = Descriptive.max(xList);
        double thresholdForColor = 0.05;
        //sometimes there is only one point
        if (vMin - vMax == 0) {
            vMin = vMax - 1;
        }
        if (hMin - hMax == 0) {
            hMin = hMax - 1;
        }
        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
                dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
        int maxYScalLen = drawAxesScale(g2d, transformer);
        //  drawAxesSciScale(g2d, transformer, true, true);
        drawAxes(g2d, title, "Relative Pyshical Position (bp)", "-Log10(p)", maxYScalLen);
        double rectangleSize = 4;
        for (int i = 0; i < pointSize; i++) {
            if (yList.getQuick(i) > (-Math.log10(thresholdForColor))) {
                g2d.setColor(Color.RED);
            } else {
                g2d.setColor(Color.BLACK);
            }
            Point2D point1 = transformer.data2ScreenPoint(xList.getQuick(i), yList.getQuick(i));
            Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                    point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
            g2d.draw(f);
            g2d.fill(f);
        }
        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);
    }

    public void drawMultiScatterPlot(List<DoubleArrayList> xList, List<DoubleArrayList> yList, List<String> legends, String title, String outputPath) throws Exception {
        if (xList == null || xList.isEmpty()) {
            System.err.println("Null  xList");
            return;
        }
        dataPlottingOffsetLeft = 40;
        int listNum = xList.size();
        for (int i = listNum - 1; i >= 0; i--) {
            if (xList.get(i) == null || xList.get(i).size() == 0) {
                System.err.println("Null xList");
                xList.remove(i);
            }
        }

        if (title != null && title.length() > 0) {
            calculateDataPlottingArea(true);
        } else {
            calculateDataPlottingArea(false);
        }

        BufferedImage image = new BufferedImage(this.canvasWidth, this.canvasHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintCanvas(g2d);

        g2d.setStroke(new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));

        //max and min value in vertical
        /*
        double vMax = Descriptive.max(yList);
        double vMin = Descriptive.min(yList);
        double hMin = Descriptive.min(xList);
        double hMax = Descriptive.max(xList);
         * 
         */
        double vMax = 1;
        double vMin = 0;
        double hMin = 0;
        double hMax = 1;

        //sometimes there is only one point
        if (vMin - vMax == 0) {
            vMin = vMax - 1;
        }
        if (hMin - hMax == 0) {
            hMin = hMax - 1;
        }
        int listSize = xList.size();
        int maxLegentWidth = g2d.getFontMetrics().stringWidth(legends.get(listSize - 1));
        for (int i = listSize - 2; i >= 0; i--) {
            int strLen = g2d.getFontMetrics().stringWidth(legends.get(i));
            if (maxLegentWidth < strLen) {
                maxLegentWidth = strLen;
            }
        }

        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(dataPlottingArea.x, dataPlottingArea.x + dataPlottingArea.width,
                dataPlottingArea.y, dataPlottingArea.y + dataPlottingArea.height, hMin, hMax, vMin, vMax);
        int maxYScalLen = drawAxesScale(g2d, transformer);
        drawAxes(g2d, title, "False positive rate", "True positive rate", maxYScalLen);

        double rectangleSize = 4;
        int strHei = g2d.getFontMetrics().getHeight();
        Color[] colors = {Color.BLUE, Color.RED, Color.ORANGE, Color.MAGENTA, Color.BLACK, Color.YELLOW};
        Point2D point0 = null;
        Line2D line = new Line2D.Double();
        for (int t = 0; t < listSize; t++) {
            int pointSize = xList.get(t).size();
            g2d.setColor(colors[t]);
            for (int i = 0; i < pointSize; i++) {
                Point2D point1 = transformer.data2ScreenPoint(xList.get(t).getQuick(i), yList.get(t).getQuick(i));
                Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                        point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                g2d.draw(f);
                g2d.fill(f);
                if (i > 0) {
                    line.setLine(point0, point1);
                    g2d.draw(line);
                }
                point0 = point1;
            }

            //g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
            g2d.setStroke(new BasicStroke(1.5f));
            double xOffset = dataPlottingArea.getWidth() / 4;
            double yOffset = dataPlottingArea.getHeight() / 2;
            g2d.drawString(legends.get(t), (float) (dataPlottingArea.x + xOffset),
                    (float) (dataPlottingArea.getY() + yOffset) + strHei * t);
            Rectangle2D.Double f = new Rectangle2D.Double(dataPlottingArea.x + xOffset + maxLegentWidth + rectangleSize,
                    dataPlottingArea.getY() + yOffset + strHei * t - rectangleSize, rectangleSize, rectangleSize);
            g2d.fill(f);
        }
        g2d.dispose();

        outputPNGFile(image, outputPath);
    }
}
