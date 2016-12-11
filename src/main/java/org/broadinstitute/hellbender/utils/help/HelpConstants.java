package org.broadinstitute.hellbender.utils.help;

public final class HelpConstants {

    public static final String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";

    public static String forumPost(String post) {
        return GATK_FORUM_URL + post;
    }

    /**
     * Definition of the group names / categories of tools.
//     * The names get parsed to make supercategories in the doc index,
//     * so be careful when making big changes -- see GATKDoclet.java toMap()
     */
    //public final static String DOCS_CAT_ENGINE = "Engine Parameters (available to all tools)";
    public final static String DOCS_CAT_SPARK = "Spark Tools";
    public final static String DOCS_CAT_SPARK_PIPELINE = "Spark Pipeline";
    public final static String DOCS_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOCS_CAT_DATA = "Sequence Data Processing Tools";
    public final static String DOCS_CAT_VARDISC = "Variant Discovery Tools";
    public final static String DOCS_CAT_VAREVAL = "Variant Evaluation Tools";
    public final static String DOCS_CAT_VARMANIP = "Variant Manipulation Tools";
    public final static String DOCS_CAT_ANNOT = "Annotation Modules";
    //public final static String DOCS_CAT_RF = "Read Filters";
    //public final static String DOCS_CAT_RODCODECS = "Resource File Codecs";
    public final static String DOCS_CAT_REFUTILS = "Reference Utilities";
    //public final static String DOCS_CAT_USRERR = "User Exceptions (Exclude)";
    public final static String DOCS_CAT_TOY = "Toy Examples (Exclude)";

    public final static String DOCS_CAT_READFILTERS = "Read Filters";


}