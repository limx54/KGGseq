/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Iterator;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class BasicPainter {
    //unit pixel

    protected int canvasWidth = 650;
    protected int canvasHeight = 400;
    protected int dataPlottingOffsetLeft = 50;
    protected int dataPlottingOffsetTop = 10;
    protected int dataPlottingOffsetBottom = 40;
    protected int dataPlottingOffsetRight = 20;
    //area to present the data point
    protected Rectangle dataPlottingArea;
    protected Color canvasBackgroundColor = Color.WHITE;
    protected Color axesColor = Color.DARK_GRAY;
    protected Color nonplotBackgroundColor = Color.LIGHT_GRAY;
    protected Font titleFont = new Font("SansSerif", Font.BOLD, 16);
    protected Font numFont = new Font("SansSerif", Font.BOLD, 12);
    protected int titleLineHeight = 30;
    protected int yScalDecimal = 2;

    public int getDataPlottingOffsetBottom() {
        return dataPlottingOffsetBottom;
    }

    public void setDataPlottingOffsetBottom(int dataPlottingOffsetBottom) {
        this.dataPlottingOffsetBottom = dataPlottingOffsetBottom;
    }

    public int getDataPlottingOffsetLeft() {
        return dataPlottingOffsetLeft;
    }

    public void setDataPlottingOffsetLeft(int dataPlottingOffsetLeft) {
        this.dataPlottingOffsetLeft = dataPlottingOffsetLeft;
    }

    public int getDataPlottingOffsetRight() {
        return dataPlottingOffsetRight;
    }

    public void setDataPlottingOffsetRight(int dataPlottingOffsetRight) {
        this.dataPlottingOffsetRight = dataPlottingOffsetRight;
    }

    public int getDataPlottingOffsetTop() {
        return dataPlottingOffsetTop;
    }

    public void setDataPlottingOffsetTop(int dataPlottingOffsetTop) {
        this.dataPlottingOffsetTop = dataPlottingOffsetTop;
    }

    public BasicPainter(int width, int height) {
        canvasWidth = width;
        canvasHeight = height;
    }

    protected void plotTriangle(Graphics2D g, int positionX, int postionY, int triangleLen) {
        Point p1 = new Point(positionX, postionY);
        Point p2 = new Point(positionX - triangleLen / 2, postionY + triangleLen);
        Point p3 = new Point(positionX + triangleLen / 2, postionY + triangleLen);

        Polygon filledPolygon = new Polygon();
        filledPolygon.addPoint(p1.x, p1.y);
        filledPolygon.addPoint(p2.x, p2.y);
        filledPolygon.addPoint(p3.x, p3.y);

        g.fillPolygon(filledPolygon);
    }

    protected void calculateDataPlottingArea(boolean drawTitle , int legendBoxWidth) throws Exception {
        if (drawTitle) {
            dataPlottingArea = new Rectangle(dataPlottingOffsetLeft, dataPlottingOffsetTop,
                    canvasWidth - dataPlottingOffsetRight - dataPlottingOffsetLeft-legendBoxWidth, canvasHeight - dataPlottingOffsetBottom - dataPlottingOffsetTop - titleLineHeight);
        } else {
            dataPlottingArea = new Rectangle(dataPlottingOffsetLeft, dataPlottingOffsetTop,
                    canvasWidth - dataPlottingOffsetRight - dataPlottingOffsetLeft-legendBoxWidth, canvasHeight - dataPlottingOffsetBottom - dataPlottingOffsetTop);
        }
    }
    
    protected void calculateDataPlottingArea(boolean drawTitle) throws Exception {
        if (drawTitle) {
            dataPlottingArea = new Rectangle(dataPlottingOffsetLeft, dataPlottingOffsetTop,
                    canvasWidth - dataPlottingOffsetRight - dataPlottingOffsetLeft, canvasHeight - dataPlottingOffsetBottom - dataPlottingOffsetTop - titleLineHeight);
        } else {
            dataPlottingArea = new Rectangle(dataPlottingOffsetLeft, dataPlottingOffsetTop,
                    canvasWidth - dataPlottingOffsetRight - dataPlottingOffsetLeft, canvasHeight - dataPlottingOffsetBottom - dataPlottingOffsetTop);
        }
    }

    /**
     *
     * @param image
     * @param outputPath
     * @throws Exception
     */
    public static void outputJPEGFile(BufferedImage image, String outputPath) throws Exception {
        ImageWriter writer = null;
        ImageTypeSpecifier type = ImageTypeSpecifier.createFromRenderedImage(image);
        Iterator iter = ImageIO.getImageWriters(type, "jpg");
        if (iter.hasNext()) {
            writer = (ImageWriter) iter.next();
        }
        if (writer == null) {
            return;
        }
        IIOImage iioImage = new IIOImage(image, null, null);
        ImageWriteParam param = writer.getDefaultWriteParam();

        param.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
        param.setCompressionQuality((float) 1.0);
        ImageOutputStream outputStream = ImageIO.createImageOutputStream(new File(outputPath));
        writer.setOutput(outputStream);
        writer.write(null, iioImage, param);
        outputStream.flush();
        outputStream.close();
        writer.dispose();
    }

    /**
     *
     * @param image
     * @param outputPath
     * @throws Exception
     */
    public static void outputPNGFile(BufferedImage image, String outputPath) throws Exception {
        ImageWriter writer = null;
        ImageTypeSpecifier type = ImageTypeSpecifier.createFromRenderedImage(image);
        Iterator iter = ImageIO.getImageWriters(type, "png");
        if (iter.hasNext()) {
            writer = (ImageWriter) iter.next();
        }
        if (writer == null) {
            return;
        }
        IIOImage iioImage = new IIOImage(image, null, null);
        ImageWriteParam param = writer.getDefaultWriteParam();
        ImageOutputStream outputStream = ImageIO.createImageOutputStream(new File(outputPath));
        writer.setOutput(outputStream);
        writer.write(null, iioImage, param);
        outputStream.flush();
        outputStream.close();
        writer.dispose();
    }

    public void paintCanvas(Graphics2D g2d) {
        g2d.clipRect(0, 0, canvasWidth, canvasHeight);
        g2d.setColor(canvasBackgroundColor);
        Rectangle2D rc = new Rectangle2D.Double((double) 0, (double) 0, (double) canvasWidth, (double) canvasHeight);
        g2d.fill(rc);
        GradientPaint pat = new GradientPaint(0f, 0f, Color.WHITE, 100f, 45f, Color.GREEN);
        g2d.setPaint(pat);

    }

    /**
    
     * @param g2
     */
    protected void drawAxes(Graphics2D g2, String title, String XLabel, String yLabel, int maxYScalLen) {
        g2.setColor(axesColor.darker());
        //vertical axis
        Line2D.Double zz = new Line2D.Double((double) dataPlottingArea.getX(),
                (double) (dataPlottingArea.getY()), (double) dataPlottingArea.getX(), (double) (dataPlottingArea.getY() + dataPlottingArea.getHeight()));

        g2.draw(zz);
        Font oldFront = g2.getFont();
        g2.setFont(numFont);
        StringBuilder sb=new StringBuilder();
        for (int i=0; i<maxYScalLen;i++) sb.append('0');
        int strYScalLen = g2.getFontMetrics().stringWidth(sb.toString());
        int strLen = g2.getFontMetrics().stringWidth(yLabel);
        int strHei = g2.getFontMetrics().getAscent();
        if (strLen > dataPlottingOffsetLeft) {
            strLen = dataPlottingOffsetLeft;
        }

        g2.rotate(-90.0 * Math.PI / 180.0);
        g2.drawString(yLabel, -(dataPlottingArea.height + strLen) / 2 - dataPlottingArea.y, dataPlottingArea.x - strYScalLen);
        g2.rotate(90.0 * Math.PI / 180.0);

        //horizontal axis
        zz = new Line2D.Double((double) dataPlottingArea.getX(),
                (double) (dataPlottingArea.getY() + dataPlottingArea.getHeight()),
                (double) (dataPlottingArea.getX() + dataPlottingArea.getWidth()),
                (double) (dataPlottingArea.getY() + dataPlottingArea.getHeight()));
        g2.draw(zz);


        strLen = g2.getFontMetrics().stringWidth(XLabel);
        g2.drawString(XLabel, dataPlottingArea.x + dataPlottingArea.width / 2 - strLen / 2,
                (float) (dataPlottingArea.getY() + dataPlottingArea.getHeight() + dataPlottingOffsetBottom * 2 / 3));

        if (title != null && title.length() > 0) {
            drawTitle(g2, title);
        }
        g2.setFont(oldFront);
    }

    protected void drawTitle(Graphics2D g2d, String title) {
        g2d.setColor(axesColor.darker());
        Font oldFront = g2d.getFont();
        if (title != null && title.length() > 0) {
            g2d.setFont(titleFont);
            int strLen = g2d.getFontMetrics().stringWidth(title);
            // ?????????
            g2d.drawString(
                    title,
                    dataPlottingArea.x + dataPlottingArea.width / 2 - strLen
                    / 2,
                    (float) (dataPlottingArea.getY()
                    + dataPlottingArea.getHeight()
                    + dataPlottingOffsetBottom + titleLineHeight / 2));
        }

        g2d.setFont(oldFront);
    }

    protected int drawAxesSciScale(Graphics2D g2, CoordinateTransformer transformer, boolean xSci, boolean ySci) {
        g2.setColor(axesColor.darker());
        String str;
        int scaleLen = 3;
        int verticalScale = 10;
        double vMin = transformer.getDataVerticalMin();
        double vMax = transformer.getDataVerticalMax();
        double hMin = transformer.getDataHorizontalMin();
        double hMax = transformer.getDataHorizontalMax();
        double gridLen = (vMax - vMin) / verticalScale;
        Font oldFront = g2.getFont();
        g2.setFont(numFont);
        int maxYScaleLen = 0;
        int str_hei = g2.getFontMetrics().getAscent();
        for (int i = 1; i <= verticalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin, vMin + i * gridLen);
            Point2D point2 = new Point2D.Double(point1.getX() + scaleLen, point1.getY());
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            if (ySci) {
                str = Util.formatPValue(vMin + i * gridLen);
            } else {
                str = Util.doubleToString(vMin + i * gridLen, 5, yScalDecimal);
            }
            if (maxYScaleLen < str.length()) {
                maxYScaleLen = str.length();
            }
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len,
                    (float) (point1.getY() + str_hei / 2));
        }

        int horizontalScale = 10;
        gridLen = (hMax - hMin) / horizontalScale;
        for (int i = 1; i <= horizontalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin + i * gridLen, vMin);
            Point2D point2 = new Point2D.Double(point1.getX(), point1.getY() - scaleLen);
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            if (xSci) {
                str = Util.formatPValue(hMin + i * gridLen);
            } else {
                str = Util.doubleToString(hMin + i * gridLen, 5, 2);
            }
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len / 2,
                    (float) (point1.getY() + str_hei));
        }
        g2.setFont(oldFront);
        return maxYScaleLen;
    }

    protected int drawAxesScale(Graphics2D g2, CoordinateTransformer transformer) {
        g2.setColor(axesColor.darker());
        String str;
        int scaleLen = 3;
        int verticalScale = 10;
        double vMin = transformer.getDataVerticalMin();
        double vMax = transformer.getDataVerticalMax();
        double hMin = transformer.getDataHorizontalMin();
        double hMax = transformer.getDataHorizontalMax();
        double gridLen = (vMax - vMin) / verticalScale;

        Font oldFront = g2.getFont();
        g2.setFont(numFont);
        int maxYScaleLen = 0;
        int str_hei = g2.getFontMetrics().getAscent();
        for (int i = 1; i <= verticalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin, vMin + i * gridLen);
            Point2D point2 = new Point2D.Double(point1.getX() + scaleLen, point1.getY());
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.doubleToString(vMin + i * gridLen, 5, yScalDecimal);
            if (maxYScaleLen < str.length()) {
                maxYScaleLen = str.length();
            }
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len,
                    (float) (point1.getY() + str_hei / 2));
        }

        int horizontalScale = 10;
        gridLen = (hMax - hMin) / horizontalScale;
        for (int i = 1; i <= horizontalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin + i * gridLen, vMin);
            Point2D point2 = new Point2D.Double(point1.getX(), point1.getY() - scaleLen);
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.doubleToString(hMin + i * gridLen, 5, 2);
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len*1/2 ,
                    (float) (point1.getY() + str_hei));
        }
        g2.setFont(oldFront);
        return maxYScaleLen;
    }

    /**
     *
     */
    public class CoordinateTransformer {

        protected int screenHorizontalMin, screenVerticalMin, screenHorizontalMax, screenVerticalMax;
        protected double dataHorizontalMin, dataVerticalMin, dataHorizontalMax, dataVerticalMax;
        double horizontalResolution, verticalResolution;

        /**
         *
         * @return
         */
        public double getDataHorizontalMax() {
            return dataHorizontalMax;
        }

        /**
         *
         * @param dataHorizontalMax
         */
        public void setDataHorizontalMax(double dataHorizontalMax) {
            this.dataHorizontalMax = dataHorizontalMax;
        }

        /**
         *
         * @return
         */
        public double getDataHorizontalMin() {
            return dataHorizontalMin;
        }

        /**
         *
         * @param dataHorizontalMin
         */
        public void setDataHorizontalMin(double dataHorizontalMin) {
            this.dataHorizontalMin = dataHorizontalMin;
        }

        /**
         *
         * @return
         */
        public double getDataVerticalMax() {
            return dataVerticalMax;
        }

        /**
         *
         * @param dataVerticalMax
         */
        public void setDataVerticalMax(double dataVerticalMax) {
            this.dataVerticalMax = dataVerticalMax;
        }

        /**
         *
         * @return
         */
        public double getDataVerticalMin() {
            return dataVerticalMin;
        }

        /**
         *
         * @param dataVerticalMin
         */
        public void setDataVerticalMin(double dataVerticalMin) {
            this.dataVerticalMin = dataVerticalMin;
        }

        /**
         *
         * @return
         */
        public double getHorizontalResolution() {
            return horizontalResolution;
        }

        /**
         *
         * @param horizontalResolution
         */
        public void setHorizontalResolution(double horizontalResolution) {
            this.horizontalResolution = horizontalResolution;
        }

        /**
         *
         * @return
         */
        public int getScreenHorizontalMax() {
            return screenHorizontalMax;
        }

        /**
         *
         * @param screenHorizontalMax
         */
        public void setScreenHorizontalMax(int screenHorizontalMax) {
            this.screenHorizontalMax = screenHorizontalMax;
        }

        /**
         *
         * @return
         */
        public int getScreenHorizontalMin() {
            return screenHorizontalMin;
        }

        /**
         *
         * @param screenHorizontalMin
         */
        public void setScreenHorizontalMin(int screenHorizontalMin) {
            this.screenHorizontalMin = screenHorizontalMin;
        }

        /**
         *
         * @return
         */
        public int getScreenVerticalMax() {
            return screenVerticalMax;
        }

        /**
         *
         * @param screenVerticalMax
         */
        public void setScreenVerticalMax(int screenVerticalMax) {
            this.screenVerticalMax = screenVerticalMax;
        }

        /**
         *
         * @return
         */
        public int getScreenVerticalMin() {
            return screenVerticalMin;
        }

        /**
         *
         * @param screenVerticalMin
         */
        public void setScreenVerticalMin(int screenVerticalMin) {
            this.screenVerticalMin = screenVerticalMin;
        }

        /**
         *
         * @return
         */
        public double getVerticalResolution() {
            return verticalResolution;
        }

        /**
         *
         * @param verticalResolution
         */
        public void setVerticalResolution(double verticalResolution) {
            this.verticalResolution = verticalResolution;
        }

        /**
         *
         * @param screenHorizontalMin
         * @param screenHorizontalMax
         * @param screenVerticalMin
         * @param screenVerticalMax
         * @param dataHorizontalMin
         * @param dataHorizontalMax
         * @param dataVerticalMin
         * @param dataVerticalMax
         */
        public void setupBasicScope(int screenHorizontalMin, int screenHorizontalMax, int screenVerticalMin, int screenVerticalMax,
                double dataHorizontalMin, double dataHorizontalMax, double dataVerticalMin, double dataVerticalMax) {
            this.screenHorizontalMin = screenHorizontalMin;
            this.screenHorizontalMax = screenHorizontalMax;
            this.screenVerticalMin = screenVerticalMin;
            this.screenVerticalMax = screenVerticalMax;

            this.dataHorizontalMin = dataHorizontalMin;
            this.dataHorizontalMax = dataHorizontalMax;
            this.dataVerticalMin = dataVerticalMin;
            this.dataVerticalMax = dataVerticalMax;

            this.horizontalResolution = (dataHorizontalMax - dataHorizontalMin) / (screenHorizontalMax - screenHorizontalMin);
            this.verticalResolution = (dataVerticalMax - dataVerticalMin) / (screenVerticalMax - screenVerticalMin);

        }
// note the vertical coordinates has been converted as we see usually

        /**
         *
         * @param Xa
         * @param Ya
         * @return
         */
        public double[] data2Screen(double Xa, double Ya) {
            double[] myout = new double[2];
            myout[0] = ((Xa - dataHorizontalMin) / horizontalResolution + screenHorizontalMin);
            myout[1] = (screenVerticalMax - (Ya - dataVerticalMin) / verticalResolution);
            return myout;
        }

        /**
         *
         * @param Xa
         * @param Ya
         * @return
         */
        public Point data2ScreenPoint(double Xa, double Ya) {
            Point myout = new Point();
            myout.setLocation(((Xa - dataHorizontalMin) / horizontalResolution + screenHorizontalMin),
                    (screenVerticalMax - (Ya - dataVerticalMin) / verticalResolution));
            return myout;
        }

        public void data2ScreenPoint(Point myout, double Xa, double Ya) {
            myout.setLocation(((Xa - dataHorizontalMin) / horizontalResolution + screenHorizontalMin),
                    (screenVerticalMax - (Ya - dataVerticalMin) / verticalResolution));
        }

        /**
         *
         * @param X6
         * @param Y6
         * @return
         */
        public double[] screen2Data(int X6, int Y6) {
            double[] myout = new double[2];
            myout[0] = (X6 - screenHorizontalMin) * horizontalResolution + dataHorizontalMin;
            myout[1] = (screenVerticalMax - Y6) * verticalResolution + dataVerticalMin;
            return myout;
        }

        /**
         *
         * @param len
         * @return
         */
        public double horizontalSegmentData2Screen(double len) {
            return len / horizontalResolution;
        }

        /**
         *
         * @param len
         * @return
         */
        public double verticalSegmentData2Screen(double len) {
            return len / verticalResolution;
        }
    }
}
