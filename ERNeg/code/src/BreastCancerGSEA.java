import java.util.*;
import org.apache.commons.io.FileUtils;
import java.io.*;
import java.util.Map.Entry;

/**
 * Basic methods used to generate status files, process output from gsea Python script and generate supplementary output tables.
 * Author: Shirley Hui
 * Date: Feb 1, 2017
 * Time: 11:26:26 AM
 */
public class BreastCancerGSEA {

    private static final  String ROOT_DIR = "/home/shirleyhui/Work/BreastCancer_PathwayAnalysis/ERNeg/";

    public  BreastCancerGSEA()
    {

    }
    public static void main(String[] args)
    {
        BreastCancerGSEA b = new BreastCancerGSEA();

        //icogs
        //b.makeSampleIxStatusFile_icogs();
        //b.makeStatusFile_icogs();

        //onco
        //b.makeStatusFile_onco();
        //b.makeSampleIxStatusFile_onco();
        
        // Run gsea python script here to create enrichment score output file
        //b.filterPathways(0.4);
 
        // Make supplementary tables
        String clusterFile = ROOT_DIR + "data/cimbabcacerneg-HISTONES-chr23-clusterbytheme.txt";
        String pathwayOutputFile = ROOT_DIR + "data/cimbabcacerneg-HISTONES-chr23-pathways.txt";
        String pesFile = ROOT_DIR + "data/pes/brca1.bcacerneg.ma.imp.all.breast_1-HISTONES-chr23.txt.pEs_Jan192016_gmt-0.4.txt";
        String emFile = ROOT_DIR + "data/em/brca1.bcacerneg.ma.imp.all.breast_1-HISTONES-chr23.txt.pEs_Jan192016_gmt-0.4.er.txt";
        b.tableS12(clusterFile, pathwayOutputFile, pesFile, emFile);
        b.tableS13(clusterFile, pathwayOutputFile, pesFile);
	b.tableS14();
    }
   public HashMap loadPathwayGenes(String pathwayToGeneFile)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(pathwayToGeneFile)));
            String line = "";
            String pw = "";
            String geneString = "";
            List geneList = new ArrayList();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!geneList.isEmpty())
                    {
                        geneString = (String)geneList.get(0);
                        for (int i = 1; i < geneList.size();i++)
                        {
                            geneString = geneString + ","+(String)geneList.get(i);

                        }
                        map.put(pw,geneString);
                    }
                    String pw0 = line.substring(1,line.length());
                    String[] splitpw = pw0.split("\t");
                    pw = splitpw[0];
                    geneList = new ArrayList();
                    geneString = "";
                }
               else
                {
                    String[] splitLine = line.split("\t");
                    String gene = splitLine[0];
                    if (!geneList.contains(gene))
                        geneList.add(gene);
                }

            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
        return map;
    }

   public HashMap loadPathwayToThemeMap(String clusterByThemeMapFilename)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(clusterByThemeMapFilename)));
            String line = "";
            String theme = "";
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    theme = line.substring(1,line.length());
                }
                else
                {
                    map.put(line,theme);
                    //System.out.println(line+"->"+theme);
                }
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
        return map;
    }


    public void tableS14()
    {

        HashMap pwToGeneMap = loadPathwayGenes(ROOT_DIR + "data/cimbabcacerneg-HISTONES-chr23-pathways-allgenes.txt");
        String clusterByThemeFile2 = ROOT_DIR + "data/cimbabcacerneg-HISTONES-chr23-clusterbytheme.txt";

        String pesFile2 = ROOT_DIR + "data/pes/brca1.bcacerneg.ma.imp.all.breast_1-HISTONES-chr23.txt.pEs_Jan192016_gmt-1.0.txt";
        String pesFile1 = ROOT_DIR + "data/pes/all.meta.gwas.bcfr.icogs.oncoarray.erpos.updMarch20161-HISTONES-chr23.txt.pEs_Jan192016_gmt-1.0.txt";

        HashMap pathwayToThemeMap2 = loadPathwayToThemeMap(clusterByThemeFile2);

        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(pesFile1)));
            String line = "";
            List pwList1 = new ArrayList();
            HashMap pwInfoMap1 = new HashMap();
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pw1 = splitLine[0];
                String es1 = splitLine[1];
                String size1 = splitLine[2];
                pwInfoMap1.put(pw1, es1 + "\t" + size1);
                pwList1.add(pw1);
            }
            System.out.println(pwList1.size());
            br.close();
            br = new BufferedReader(new FileReader(new File(pesFile2)));
            line = "";
            List pwList2 = new ArrayList();
            HashMap pwInfoMap2 = new HashMap();

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pw2 = splitLine[0];
                String es2 = splitLine[1];
                String size2 = splitLine[2];
                pwList2.add(pw2);
                pwInfoMap2.put(pw2,es2 + "\t" + size2);

            }
            br.close();
            String outfilename = "TableS14.txt";
            BufferedWriter bw =  new BufferedWriter(new FileWriter(new File(outfilename)));

            for (int i = 0 ; i < pwList2.size();i++)
            {
                String pw2 = (String)pwList2.get(i);
                String pwInfo2 = (String)pwInfoMap2.get(pw2);
                String[] splitInfo2 = pwInfo2.split("\t");
                String es2 = splitInfo2[0];
                String size2 = splitInfo2[1];
                if (!pwList1.contains(pw2))
                {
                    String theme = (String)pathwayToThemeMap2.get(pw2);
                    String geneString = (String) pwToGeneMap.get(pw2);
                    if (theme== null)
                    {
                        continue;
                    }
                    String[] splitpw = pw2.split("%");
                    System.out.println(theme + "\t" + splitpw[0] +"\t"+ splitpw[1] + "\t" + splitpw[2] + "\t" + es2+ "\t" + size2+"\t" + geneString) ;
                    bw.write(theme+ "\t" + splitpw[0] +"\t"+ splitpw[1] + "\t" + splitpw[2] + "\t" + es2 + "\t" + size2 +"\t" + geneString +"\n") ;

                }
            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
    }

    public void tableS12(String clusterFile, String pathwayOutputFile, String pesFile, String emFile)
    {
        try
        {
            HashMap finalMap = new HashMap();
            HashMap map = new HashMap();
            BufferedReader br = new BufferedReader(new FileReader(new File(emFile)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                map.put(splitLine[0],splitLine[2]);
            }
            br.close();

            HashMap clusterToPathways = new HashMap();
            br = new BufferedReader(new FileReader(new File(clusterFile)));
            line = "";
            List pathways = new ArrayList();
            String clusterName = "";
            HashMap pathwaysToCluster= new HashMap();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!pathways.isEmpty())
                    {
                        clusterToPathways.put(clusterName,pathways);
                    }
                    String[] splitLine = line.split("\t");

                    clusterName = splitLine[0].substring(1,splitLine[0].length());
                    pathways = new ArrayList();
                }
                else
                {
                    pathways.add(line);
                    pathwaysToCluster.put(line,clusterName);
                }
            }
            clusterToPathways.put(clusterName,pathways);

            br.close();

            HashMap pathwayToPESSizeMap = new HashMap();
            br = new BufferedReader(new FileReader(new File(pesFile)));
            line = "";

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathwayName = splitLine[0];
                String pes = splitLine[1];
                pathwayToPESSizeMap.put(pathwayName, pes);
            }
            br.close();
            br = new BufferedReader(new FileReader(new File(pathwayOutputFile)));
            line = "";
            String pathwayName = "";
            List genes = new ArrayList();
            clusterName = "";
            HashMap pathwayToGenes = new HashMap();
            HashMap pathwayToFullName = new HashMap();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!genes.isEmpty())
                    {
                        pathwayToGenes.put(pathwayName, genes);
                        genes = new ArrayList();
                    }
                    String[] splitLine =line.split("\t");
                    pathwayName = splitLine[0].substring(1,splitLine[0].length());
                    pathwayToFullName.put(pathwayName, line.substring(1,line.length())) ;
                }
                else
                {
                    genes.add(line);
                }
            }
            br.close();
            Set keys = clusterToPathways.keySet();
            List keyList = new ArrayList(keys);
            for (int i = 0; i < keyList.size();i++)
            {
                String clusterName0 = (String) keyList.get(i);
                List pathways0 = (List) clusterToPathways.get(clusterName0);

                for(int j = 0;j < pathways0.size();j++)
                {
                    String pathway0 = (String)pathways0.get(j);
                    String pathway00 = (String)pathwayToFullName.get(pathway0);
                    if (pathway00 == null)
                        continue;
                    String pes00 = (String)pathwayToPESSizeMap.get(pathway0);
                    if (pes00 == null)
                        continue;
                    String[] splitpathway00 = pathway00.split("\t");

                    //System.out.println("#Genes: " +splitpathway00[1]);
                    //System.out.println("Score (ES): " +  pes00);
                    //System.out.println("Gene\t#SNPS\tMost Sig. SNP\tP-value\tDistance to Gene (bp)");
                    String ses = (String)map.get(pathway0);
                    List genes0 = (List)pathwayToGenes.get(pathway0);
                    for (int k = 0 ; k < genes0.size();k++)
                    {
                        //System.out.println(genes0.get(k));
                        String pw = splitpathway00[0];
                        String[] splitPw = pw.split("%");
                        //System.out.println(clusterName0 + "\t" +splitPw[0] + "\t" + splitPw[1] + "\t" + splitPw[2] + "\t" + pes00 + "\t"+ ses + "\t"+ genes0.get(k));
                        List info = (List)finalMap.get(clusterName0);
                        if (info == null)
                           info = new ArrayList();
                        info.add(clusterName0 + "\t" +splitPw[0] + "\t" + splitPw[1] + "\t" + splitPw[2] + "\t" + pes00 + "\t"+ ses + "\t"+ genes0.get(k));
                        finalMap.put(clusterName0,info);
                    }
                }

            }
            System.out.println("######");
            keys = finalMap.keySet();
            keyList = new ArrayList(keys);
            Collections.sort(keyList);
            String outfilename = "TableS12.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            for (int i = 0; i < keyList.size();i++)
            {
                String k = (String)keyList.get(i);
                List info = (List)finalMap.get(k);
                for (int j =0;j < info.size();j++)
                {
                    //System.out.println(info.get(j));
                    bw.write(info.get(j)+"\n");
                }
            }
            bw.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
    }

    public void tableS13(String clusterFile, String pathwayOutputFile, String pesFile)
    {
        try
        {
            HashMap clusterToPathways = new HashMap();
            BufferedReader br = new BufferedReader(new FileReader(new File(clusterFile)));
            String line = "";
            List pathways = new ArrayList();
            String clusterName = "";
            HashMap pathwaysToCluster= new HashMap();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!pathways.isEmpty())
                    {
                        clusterToPathways.put(clusterName,pathways);
                    }
                    String[] splitLine = line.split("\t");
                    clusterName = splitLine[0].substring(1,splitLine[0].length());
                    pathways = new ArrayList();
                }
                else
                {
                    pathways.add(line);
                    pathwaysToCluster.put(line,clusterName);
                }
            }
            clusterToPathways.put(clusterName,pathways);

            br.close();

            HashMap pathwayToPESSizeMap = new HashMap();
            br = new BufferedReader(new FileReader(new File(pesFile)));
            line = "";

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathwayName = splitLine[0];
                String pes = splitLine[1];
                pathwayToPESSizeMap.put(pathwayName, pes);
            }

            br.close();
            br = new BufferedReader(new FileReader(new File(pathwayOutputFile)));
            line = "";
            String pathwayName = "";
            List genes = new ArrayList();
            clusterName = "";
            HashMap pathwayToGenes = new HashMap();
            HashMap pathwayToFullName = new HashMap();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!genes.isEmpty())
                    {
                        pathwayToGenes.put(pathwayName, genes);
                        genes = new ArrayList();
                    }
                    String[] splitLine =line.split("\t");
                    pathwayName = splitLine[0].substring(1,splitLine[0].length());
                    pathwayToFullName.put(pathwayName, line.substring(1,line.length())) ;
                }
                else
                {
                    genes.add(line);
                }
            }
            br.close();

            Set keys = clusterToPathways.keySet();
            List keyList = new ArrayList(keys);
            HashMap map = new HashMap();
            for (int i = 0; i < keyList.size();i++)
            {
                String clusterName0 = (String) keyList.get(i);
                //System.out.println("##################################################");
                //System.out.println("##"+clusterName0);
                //System.out.println("---------------------------------------------------");
                //System.out.println("Theme: "+clusterName0);
                List pathways0 = (List) clusterToPathways.get(clusterName0);
                List uniqueGenes = new ArrayList();
                HashMap valuesToGeneList =new HashMap();
                for(int j = 0;j < pathways0.size();j++)
                {
                    String pathway0 = (String)pathways0.get(j);
                    //System.out.println("***"+pathway0);
                    List genes0 = (List)pathwayToGenes.get(pathway0);
                    if (genes0 == null)
                        continue;
                    for (int k = 0 ; k < genes0.size();k++)
                    {
                        if (!uniqueGenes.contains(genes0.get(k)))
                        {
                            String[] splitGene = ((String)genes0.get(k)).split("\t");
                            double pvalue= Double.parseDouble(splitGene[3]);
                            List geneList = (List)valuesToGeneList.get(pvalue);
                            if (geneList == null)
                                geneList = new ArrayList();
                            if (!geneList.contains(genes0.get(k)))
                            {
                                geneList.add(genes0.get(k));
                                valuesToGeneList.put(pvalue, geneList);
                                uniqueGenes.add(genes0.get(k));
                            }

                        }
                    }
                }
                Set keys0 = valuesToGeneList.keySet();
                List keyList0 = new ArrayList(keys0);
                Collections.sort(keyList0);

                //System.out.println("##UNIQUE GENES####################################");
                //System.out.println("Unique Genes:");
                //System.out.println("Gene\t#SNPS\tMost Sig. SNP\tP-value\tDistance to Gene (bp)");
                for (int j = 0; j < keyList0.size();j++)
                {
                    double pvalue = (Double)keyList0.get(j);
                   List geneList = (List)valuesToGeneList.get(pvalue);
                    for(int k = 0; k < geneList.size();k++)
                    {
                        //System.out.println(uniqueGenes.get(j));
                        //System.out.println(clusterName0 +"\t"+geneList.get(k));
                        //bw.write(clusterName0 +"\t"+geneList.get(k) + "\n");
                        List geneList0 = (List)map.get(clusterName0);
                        if (geneList0==null)
                            geneList0 = new ArrayList();
                        geneList0.add(geneList.get(k));
                        map.put(clusterName0,geneList0);

                    }
                }
            }
            String outfilename = "TableS13.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            keys = map.keySet();
            keyList = new ArrayList(keys);
            Collections.sort(keyList);
            for (int i = 0; i < keyList.size();i++)
            {
                String cluster = (String)keyList.get(i);
                List geneList = (List)map.get(cluster);
                for(int k = 0; k < geneList.size();k++)
                {
                    System.out.println(cluster +"\t"+geneList.get(k));
                    bw.write(cluster +"\t"+geneList.get(k) + "\n");
                }
          }
          bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
    }
   
    public void tableS1(String clusterFile, String emFile)
    {
        try
        {
            HashMap map = new HashMap();
            BufferedReader br = new BufferedReader(new FileReader(new File(emFile)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                map.put(splitLine[0],splitLine[2]);
            }
            br.close();

            HashMap clusterToPathways = new HashMap();
            br = new BufferedReader(new FileReader(new File(clusterFile)));
            line = "";
            List pathways = new ArrayList();
            String clusterName = "";
            HashMap pathwaysToCluster= new HashMap();
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    if (!pathways.isEmpty())
                    {
                        clusterToPathways.put(clusterName,pathways);
                    }
                    String[] splitLine = line.split("\t");

                    clusterName = splitLine[0].substring(1,splitLine[0].length());
                    pathways = new ArrayList();
                }
                else
                {
                    pathways.add(line);

                    pathwaysToCluster.put(line,clusterName);
                }
            }
            clusterToPathways.put(clusterName,pathways);

            br.close();

            Set keys = clusterToPathways.keySet();
            List keyList = new ArrayList(keys);
            Collections.sort(keyList);
            String outfilename = "test.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            for (int i = 0; i < keyList.size();i++)
            {
                String k = (String)keyList.get(i);
                List pwList = (List)clusterToPathways.get(k);
                for (int j =0;j < pwList.size();j++)
                {
                    String pw = (String)pwList.get(j);
                    String[] splitPw = pw.split("%");

                    bw.write(k +"\t"+ splitPw[0] +"\t"+splitPw[1] + "\t"+ splitPw[2]+"\n");
                }
            }
            bw.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
    }

    public List getAshaTopSigPathways()
    {
        String filename = ROOT_DIR + "data/AshaSigPathways-erneg.txt";

        List pathways = new ArrayList();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String [] splitLine= line.split("\t");
                String pathway = splitLine[0].toUpperCase()+"%" + splitLine[4] +"%" + splitLine[5];
                pathway = pathway.toUpperCase();
                pathways.add(pathway);
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();

        }
        return pathways;
    }
    public List getInGMTPathways(List pathways)
    {
        String gmtFilename = ROOT_DIR + "Human_GOBP_AllPathways_no_GO_iea_January_19_2016_symbol.gmt";
        List gmtPathways = new ArrayList();
        List retList = new ArrayList();

        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(gmtFilename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String [] splitLine= line.split("\t");
                String pathway = splitLine[0];
                gmtPathways.add(pathway);
            }
            br.close();
            for (int i = 0; i < pathways.size();i++)
            {
                String pathway = (String)pathways.get(i);

                if (gmtPathways.contains(pathway))
                {
                    retList.add(pathway);
                }
               else
                {
                    System.out.println("*" +pathway);
                }
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();

        }
        return retList;
    }

    public void filterPathways(double cutoff)
    {
        List sortedList = new ArrayList();
        System.out.println("===== Sort pathways ===");
        System.out.println("TPR cutoff: " + cutoff);
        List ashaPathways = getAshaTopSigPathways();
        int n = ashaPathways.size();
        System.out.println("Num asha pathways: " + ashaPathways.size());
        ashaPathways = getInGMTPathways(ashaPathways);
        System.out.println("Num not in GMT: " + (n-ashaPathways.size()));
        int numAshaPathways = ashaPathways.size();
        System.out.println("Num asha pathways: " + numAshaPathways);

        int numToBreak = (int)((double)numAshaPathways*cutoff);
        System.out.println("Num to break: " + numToBreak);
        try
        {
            String filename = ROOT_DIR + "data/pes/all.meta.gwas.icogs.oncoarray.erneg.nogenomic1.txt.pEs_Jan192016_gmt-log2p.txt";
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            HashMap map = new HashMap();
            int num =0;
            int num2 =0;
            int num3 = 0;
            HashMap lineMap = new HashMap();
            HashMap pathwayCount = new HashMap();
            List pathwaysSeen = new ArrayList();
            while((line=br.readLine())!=null)
            {
                String [] splitLine= line.split("\t");
                String pathway = splitLine[0];
                int size = Integer.parseInt(splitLine[2]);
                double ES = Double.parseDouble(splitLine[1]);

                if (ES < 0)
                {
                    num3 = num3+1;
                    if (ashaPathways.contains(pathway))
                    {
                        if (!pathwaysSeen.contains(pathway))
                        {
                            num2 = num2+1;
                            pathwaysSeen.add(pathway);

                        }
                        Integer count = (Integer)pathwayCount.get(pathway);
                        if (count == null)
                            count = 0;
                        count = count +1;
                        pathwayCount.put(pathway, count);

                    }
                    continue;
                }
                if (size < 10 || size > 200)
                {
                    num= num+1;
                    if (ashaPathways.contains(pathway))
                    {
                        if (!pathwaysSeen.contains(pathway))
                        {
                            num2 = num2+1;
                            pathwaysSeen.add(pathway);

                        }
                        Integer count = (Integer)pathwayCount.get(pathway);
                        if (count == null)
                            count = 0;
                        count = count +1;
                        pathwayCount.put(pathway, count);
                    }

                    continue;
                }

                map.put(pathway,ES);
                lineMap.put(pathway, line);
            }
            br.close();

            Set keys = pathwayCount.keySet();
            List keyList = new ArrayList(keys);
            for (int i=0; i < keyList.size();i++)
            {
                String key = (String)keyList.get(i);

                int count = (Integer)pathwayCount.get(key);
                if (count > 1)
                    System.out.println(key + "=" + count);
            }

            System.out.println("Num < 10 or > 200: " + num);
            System.out.println("Num < 0: " + num3);

            System.out.println("Num asha removed:" + num2);
            System.out.println("Num in map: " + map.keySet().size());

            List<Entry<String, Double>> list = new LinkedList<Entry<String, Double>>(map.entrySet());

            // Sorting the list based on values
            Collections.sort(list, new Comparator<Entry<String, Double>>()
            {
                public int compare(Entry<String, Double> o1,
                                   Entry<String, Double> o2)
                {
                    return o2.getValue().compareTo(o1.getValue());
                }
            });

            // Maintaining insertion order with the help of LinkedList
            int hit = 0;
            for (Entry<String, Double> entry : list)
            {
                //System.out.println(entry.getKey() + "," + entry.getValue());
                String pathway =    entry.getKey();
                if (ashaPathways.contains(pathway))
                {
                    hit = hit +1;
                }

                sortedList.add(pathway);

                if (hit == numToBreak)
                {
                    System.out.println("Breaking at : " + hit);
                    break;
                }

            }

            System.out.println("Number TP pathways: "+ hit);
            System.out.println("Number FP pathways:" + (sortedList.size() - hit));

	    String outfilename = ROOT_DIR +  "data/pes/all.meta.gwas.icogs.oncoarray.erneg.nogenomic1.txt.pEs_Jan192016_gmt-"+cutoff+".txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            for (int i = 0;i < sortedList.size();i++)
            {
                String pathway = (String)sortedList.get(i);
                String line0 = (String)lineMap.get(pathway);
                bw.write(line0+"\n");

            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }
    public void makeEnrichmentFile(double cutoff, String type, String gene)
    {
        try
        {
            String filename = ROOT_DIR + "data/pes/all.meta.gwas.icogs.oncoarray.erneg.nogenomic1.txt.pEs_Jan192016_gmt-"+cutoff+".txt";

            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            double max = Double.MIN_VALUE;
            double min = Double.MAX_VALUE;
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                double ES = Double.parseDouble(splitLine[1]);
                if (ES <= 0)
                    continue;
                if (ES > max)
                    max = ES;
                if (ES < min)
                    min = ES;
            }
            br.close();

            br = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            String outfilename = ROOT_DIR + "data/em/all.meta.gwas.icogs.oncoarray.erneg.nogenomic1.txt.pEs_Jan192016_gmt-"+cutoff+".er.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            bw.write("Pathway Id\tDescription\tpvalue\t\t\n");
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                double ES = Double.parseDouble(splitLine[1]);
                if (ES <= 0)
                {
                    bw.write(pathway + "\tNA\t" + ES+"\n");

                }
                else
                {
                    String size = splitLine[2];
                    double sES = (ES-min)/(max-min);
                    bw.write(pathway + "\tNA\t" + sES+"\n");
                }
            }
            br.close();
            bw.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();

        }
    }

    public void makeSampleIxStatusFile_icogs()
    {
        List colIxList = new ArrayList();
        List oncIdList = new ArrayList();

        try
        {
            String line = "";

            String filename = ROOT_DIR + "data/pheno/B1_iCOGS_phenotype_distribution_080216-cleaned.csv";
            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            String header = br.readLine();
            String[] splitHeader = header.split(",");
            int ixx = 0;
            for (int i =0; i < splitHeader.length;i++)
            {
                String var = splitHeader[i];
                if (var.equals("censvar"))
                {
                    System.out.println("**" + ixx);
                    break;
                }
                ixx = ixx+1;
            }

            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split(",");
                int oncId = Integer.parseInt(splitLine[2]);
                oncIdList.add(oncId);
            }
            br.close();

            filename = ROOT_DIR + "data/sample_orders/sample_order_icogs_brca1.txt";
            br  = new BufferedReader(new FileReader(new File(filename)));
            int ix = 1;
            HashMap map = new HashMap();
            List sampleIdList = new ArrayList();
            while((line = br.readLine())!=null)
            {
                int oncId = Integer.parseInt(line);
                if (oncIdList.contains(oncId))
                {
                    colIxList.add(ix);
                    map.put(oncId, ix);
                    sampleIdList.add(oncId);
                }
                ix = ix+1;
            }
            System.out.println("Number of col ixs: " + colIxList.size());
            br.close();

            String outfilename = ROOT_DIR + "data/sample_orders/brca1_icogs_sample_ix_status.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            filename = ROOT_DIR + "data/pheno/B1_iCOGS_phenotype_distribution_080216-cleaned.csv";
            br  = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            br.readLine();
            int num = 0;
            while((line = br.readLine())!=null)
            {
                String [] splitLine = line.split(",");
                int oncId = Integer.parseInt(splitLine[2]);
                String status = splitLine[ixx];
                if (!status.equals("0") && !status.equals("1"))
                    System.out.println("* "+oncId + "\t" + status);
                if (sampleIdList.contains(oncId))
                {

                    int colIx = (Integer)map.get(oncId);
                    if (status.equals("0") || status.equals("1"))
                    {
                        bw.write(oncId +"\t"+ colIx+ "\t" + status + "\n");
                        num = num +1;
                    }
                    else if (status.equals("1") || status.equals("2") || status.equals("3"))
                    {
                        bw.write(oncId +"\t"+ colIx+ "\t1\n");
                        num = num +1;

                    }
                }
            }
            System.out.println("Num written: " + num);
            bw.close();
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public void makeSampleIxStatusFile_onco()
    {
        List colIxList = new ArrayList();
        List oncIdList = new ArrayList();

        try
        {
            String line = "";

            String filename = ROOT_DIR + "data/pheno/B1_Onco_phenotype_distribution_311215-cleaned.csv";

            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            String header = br.readLine();
            String[] splitHeader = header.split(",");
            int ixx = 0;
            for (int i =0; i < splitHeader.length;i++)
            {
                String var = splitHeader[i];
                if (var.equals("censvar"))
                {
                    System.out.println("**" + ixx);
                    break;
                }
                ixx = ixx+1;
            }

            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split(",");
                int oncId = Integer.parseInt(splitLine[2]);
                oncIdList.add(oncId);
            }
            br.close();


            filename = ROOT_DIR + "data/sample_orders/sample_order_onco_brca1.txt";

            br  = new BufferedReader(new FileReader(new File(filename)));
            int ix = 1;
            HashMap map = new HashMap();
            List sampleIdList = new ArrayList();
            while((line = br.readLine())!=null)
            {
                int oncId = Integer.parseInt(line);
                if (oncIdList.contains(oncId))
                {
                    colIxList.add(ix);
                    map.put(oncId, ix);
                    sampleIdList.add(oncId);
                }
                ix = ix+1;
            }
            System.out.println("Number of col ixs: " + colIxList.size());
            br.close();

            String outfilename = ROOT_DIR + "data/sample_orders/brca1_onco_sample_ix_status.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            filename = ROOT_DIR + "data/pheno/B1_Onco_phenotype_distribution_311215-cleaned.csv";
            br  = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            br.readLine();
            int num = 0;
            while((line = br.readLine())!=null)
            {
                String [] splitLine = line.split(",");
                int oncId = Integer.parseInt(splitLine[2]);
                String status = splitLine[ixx];
                if (!status.equals("0") && !status.equals("1"))
                    System.out.println("* "+oncId + "\t" + status);
                if (sampleIdList.contains(oncId))
                {

                    int colIx = (Integer)map.get(oncId);
                    if (status.equals("0") || status.equals("1"))
                    {
                        bw.write(oncId +"\t"+ colIx+ "\t" + status + "\n");
                        num = num +1;
                    }
                    else if (status.equals("1") || status.equals("2") || status.equals("3"))
                    {
                        bw.write(oncId +"\t"+ colIx+ "\t1\n");
                        num = num +1;

                    }
                }
            }
            System.out.println("Num written: " + num);
            bw.close();
            br.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }


    public void makeStatusFile_onco()
    {
        try
        {
            int maxIx =  -1;
            String filename = ROOT_DIR + "data/sample_orders/brca1_onco_sample_ix_status.txt";
            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            HashMap map = new HashMap();
            HashMap map2 = new HashMap();

            String line = "";
            int num =0;
            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String oncId = splitLine[0];
                int ix = Integer.parseInt(splitLine[1]);
                if (ix > maxIx)
                    maxIx = ix;
                map.put(ix, oncId);
                map2.put(oncId, splitLine[2]);
                num = num +1;
            }
            br.close();


            String outfilename = ROOT_DIR + "data/status/onco/Onco_euro_status0.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            for (int i = 0 ; i < maxIx+1;i++)
            {
                String oncId = (String)map.get(i);
                if (oncId == null)
                    continue;
                String status = (String)map2.get(oncId);
                if (status == null)
                    System.out.println(oncId + "\t"+status);
                bw.write(oncId +"\t"+status +"\n");
            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }
    public void makeStatusFile_icogs()
    {
        try
        {
            int maxIx =  -1;
            String filename = ROOT_DIR + "data/sample_orders/brca1_icogs_sample_ix_status.txt";
            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            HashMap map = new HashMap();
            HashMap map2 = new HashMap();

            String line = "";
            int num =0;
            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String oncId = splitLine[0];
                int ix = Integer.parseInt(splitLine[1]);
                if (ix > maxIx)
                    maxIx = ix;
                map.put(ix, oncId);
                map2.put(oncId, splitLine[2]);
                num = num +1;
            }
            br.close();


            String outfilename = ROOT_DIR + "data/status/icogs/Icogs_euro_status0.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            for (int i = 0 ; i < maxIx+1;i++)
            {
                String oncId = (String)map.get(i);
                if (oncId == null)
                    continue;
                String status = (String)map2.get(oncId);
                if (status == null)
                    System.out.println(oncId + "\t"+status);
                bw.write(oncId +"\t"+status +"\n");
            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public void removeRedundantHistoneGenesFromSnpFile(String type)
    {
        try
        {
            String snpMapFile = ROOT_DIR + "refdata/snps_to_genes_brca1.bcacerneg.ma.imp.all.breast_p5e-2.predicted.targets.txt";
            String outSnpMapFile = ROOT_DIR + "refdata/snps_to_genes_brca1.bcacerneg.ma.imp.all.breast-HISTONES-p5e-2.predicted.targets.txt";

            BufferedReader br = new BufferedReader(new FileReader(new File(snpMapFile)));
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outSnpMapFile)));

            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine= line.split("\t");
                String snp = splitLine[0];
                String gene = splitLine[1];
                String [] splitGene = gene.split(",");
                boolean hist= false;
                for (int i = 0; i < splitGene.length;i++)
                {
                    String gene1 = splitGene[i];
                    if (gene1.startsWith("HIST"))
                    {
                        hist = true;
                        break;
                    }
                }
                if (hist)
                {
                    if (type.equals("overall"))
                    {
                        if (!snp.equals("6_27657944_G_A") &&  // HIST1
                                !snp.equals("1_149219841_C_G") &&   //HIST2H2AB
                                !snp.equals("1_228517406_A_G") &&   // HIST3
                                !snp.equals("12_14871747_C_CTAT")) // HIST4
                            continue;
                    }
                    else if (type.equals("erpos"))
                    {
                        if (!snp.equals("6_26182960_T_C") &&  // HIST1
                                !snp.equals("1_120905164_G_C") &&   // HIST2H2BA
                                !snp.equals("1_149864270_G_A") &&   //HIST2H2AB
                                !snp.equals("1_228614743_T_C") &&   // HIST3
                                !snp.equals("12_14921974_T_A")) // HIST4
                            continue;

                    }
                    else if (type.equals("cimbabcacerneg"))
                    {
                        if (!snp.equals("6_27725030_G_T") &&  // HIST1H4I
                                !snp.equals("1_149876124_T_C") &&   // HIST2, NOTE include SV2A
                                !snp.equals("1_228472982_G_A") &&   // HIST3
                                !snp.equals("12_15074570_G_A")) // HIST4
                            continue;


                    }
                }
                bw.write(line+"\n");


            }
            br.close();
            bw.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception:  "+ e);
            e.printStackTrace();
        }
    }
}
