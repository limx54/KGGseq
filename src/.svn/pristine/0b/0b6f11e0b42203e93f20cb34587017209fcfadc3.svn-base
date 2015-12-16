/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.randomforests;

import java.io.Serializable;

/**
 *
 * @author mxli
 */
public class MyRandomTree implements Serializable {

    private static final long serialVersionUID = 300L;
    int classNum = 2;
     int m_Attribute = -1;
    /** The subtrees appended to this tree (node). */
    protected MyRandomTree[] m_Successors;
    /** The proportions of training instances going down each branch. */
    protected double[] m_Prop = null;
    /** Class probabilities from the training vals. */
    protected double[] m_ClassProbs = null;
    /** The split point. */
    protected double m_SplitPoint = Double.NaN;

    public int getM_Attribute() {
        return m_Attribute;
    }

    public void setM_Attribute(int m_Attribute) {
        this.m_Attribute = m_Attribute;
    }

    public double[] getM_ClassProbs() {
        return m_ClassProbs;
    }

    public void setM_ClassProbs(double[] m_ClassProbs) {
        this.m_ClassProbs = m_ClassProbs;
    }

    public double[] getM_Prop() {
        return m_Prop;
    }

    public void setM_Prop(double[] m_Prop) {
        this.m_Prop = m_Prop;
    }

    public double getM_SplitPoint() {
        return m_SplitPoint;
    }

    public void setM_SplitPoint(double m_SplitPoint) {
        this.m_SplitPoint = m_SplitPoint;
    }

    public MyRandomTree[] getM_Successors() {
        return m_Successors;
    }

    public void setM_Successors(MyRandomTree[] m_Successors) {
        this.m_Successors = m_Successors;
    }
    /**
     * Computes class distribution of an instance using the FastRandomTree.<p>
     *
     * In Weka's RandomTree, the distributions were normalized so that all
     * probabilities sum to 1; this would abolish the effect of instance weights
     * on voting. In FastRandomForest 0.97 onwards, the distributions are
     * normalized by dividing with the number of instances going into a leaf.<p>
     * 
     * @param instance the instance to compute the distribution for
     * @return the computed class distribution
     * @throws Exception if computation fails
     */
    
    public double[] distributionForInstance(double[] instance) throws Exception {
        double[] returnedDist = null;

        if (m_Attribute > -1) {  // ============================ node is not a leaf

            if (Double.isNaN(instance[m_Attribute])) {  // ---------------- missing value
                returnedDist = new double[classNum];
                // split instance up
                for (int i = 0; i < m_Successors.length; i++) {
                    double[] help = m_Successors[i].distributionForInstance(instance);
                    if (help != null) {
                        for (int j = 0; j < help.length; j++) {
                            returnedDist[j] += m_Prop[i] * help[j];
                        }
                    }
                }

            } else { // ------------------------------------------ numeric attributes
                if (instance[m_Attribute] <= m_SplitPoint) {
                    returnedDist = m_Successors[0].distributionForInstance(instance);
                } else {
                    returnedDist = m_Successors[1].distributionForInstance(instance);
                }
            }
            return returnedDist;
        } else { // =============================================== node is a leaf
            return m_ClassProbs;
        }
    }
}
