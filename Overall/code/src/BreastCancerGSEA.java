import java.util.*;
import org.apache.commons.io.FileUtils;
import java.io.*;

/**
 * Utility functions for Breast Cancer data
 * Author: Shirley Hui
 * Date: Feb 1, 2017
 * Time: 11:26:26 AM
 */
public class BreastCancerGSEA {

    private static final String DATA_DIR = "/home/shirleyhui/Work/BreastCancer_PathwayAnalysis/Overall/data/";
    //private static final  String ROOT_DIR = "/home/shirleyhui/Work/BreastCancer_PathwayAnalysis/Overall/";

    public  BreastCancerGSEA()
    {

    }

    public static void main(String[] args)
    {
        BreastCancerGSEA b = new BreastCancerGSEA();
        //onco
        //b.makeSampleIxStatusFile_onco();
        //b.makeStatusFile_onco();
        //b.stitchPermOutFiles("onco","0.05");

        //icogs
        //b.makeSampleIxStatusFile_icogs();
        //b.makeStatusFile_icogs();
        //b.stitchPermOutFiles("icogs","0.05");

        // Run python scripts to compute enrichment scores and FDRs before running below

        //b.makeEnrichmentFileFDRFilterSingleSigGene(0.05);

        String clusterFilename = "overall-predtargets-pathways-clusterall-fdr-0.05-or.txt";
        String pesFilename = "all.meta.gwas.icogs.oncoarray.overall.nogenomiccontrol1-HISTONES-5e-2.predicted.targets.txt.pEs_Jan192016_gmt-log2p-qcsnpsexcluded.txt";
        String pathwayOutputFilename = "overall-predtargets-pathways-siggenes-qcexcluded.txt";
        //b.tableS22(clusterFilename, pathwayOutputFilename, pesFilename);
        b.tableS23(clusterFilename, pathwayOutputFilename, pesFilename);
        //b.writePesFDRFile_or(pesFilename, clusterFilename, 0.05);
        //b.writeSigPathwaysFDR_or(0.05);

    }

    public void makeSampleIxStatusFile_icogs()
    {
        List colIxList = new ArrayList();
        List oncIdList = new ArrayList();

        try
        {
            String line = "";

            String filename = DATA_DIR + "pheno/iCOGS_all_euro_pheno_v2_bader_345.txt";
            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            br.readLine();
            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                int oncId = Integer.parseInt(splitLine[0]);
                oncIdList.add(oncId);
            }
            br.close();

            filename = DATA_DIR + "sample_orders/sample_order_bcac_icogs.txt";
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

            String outfilename = DATA_DIR + "sample_orders/bcac_icogs_sample_ix_status.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            filename = DATA_DIR + "pheno/iCOGS_all_euro_pheno_v2_bader_345.txt";
            br  = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            br.readLine();
            int num = 0;
            while((line = br.readLine())!=null)
            {
                String [] splitLine = line.split("\t");
                int oncId = Integer.parseInt(splitLine[0]);
                String status = splitLine[4];
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

            String filename = DATA_DIR + "pheno/Onco_euro_excluding_iCOGS_euro_overlap_pheno_v6_bader_345.txt";
            BufferedReader br  = new BufferedReader(new FileReader(new File(filename)));
            br.readLine();
            while((line = br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                int oncId = Integer.parseInt(splitLine[0]);
                oncIdList.add(oncId);
            }
            br.close();
            filename = DATA_DIR + "sample_orders/sample_order_bcac_onco.txt";

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

            String outfilename = DATA_DIR + "sample_orders/bcac_onco_sample_ix_status.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            filename = DATA_DIR + "pheno/Onco_euro_excluding_iCOGS_euro_overlap_pheno_v6_bader_345.txt";
            br  = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            br.readLine();
            int num = 0;
            while((line = br.readLine())!=null)
            {
                String [] splitLine = line.split("\t");
                int oncId = Integer.parseInt(splitLine[0]);
                String status = splitLine[5];
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
            String filename = DATA_DIR + "sample_orders/bcac_onco_sample_ix_status.txt";
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


            String outfilename = DATA_DIR + "status/onco/Onco_euro_status0.txt";
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
            String filename = DATA_DIR + "sample_orders/bcac_icogs_sample_ix_status.txt";
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


            String outfilename = DATA_DIR + "status/icogs/Icogs_euro_status0.txt";
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

    public void stitchPermOutFiles(String type, String range)
    {
        try
        {   System.out.println("Stitching chisq out files into perm files...: " + type + "," + range);
            List brList = new ArrayList();
            List bwList = new ArrayList();
            File dirFile = new File(DATA_DIR + "bcac_"+type+"_chisq/"+range+"/out/");
            File[] files = dirFile.listFiles();

            for (int i = 0 ; i < files.length;i++)
            {
                BufferedReader br = new BufferedReader(new FileReader(files[i]));
                brList.add(br);
            }

            for (int i =0;i<1001;i++)
            {
                String outfilename = DATA_DIR +"bcac_"+type+"_chisq/"+range+"/perm/bcac_"+type+"_euro_dosages_"+range+".perm"+i+".txt" ;
                bwList.add(outfilename);
            }
            System.out.println("Num outfiles added: " + bwList.size());
            for (int i = 0; i < brList.size();i++)
            {
                BufferedReader br = (BufferedReader)brList.get(i);
                String line = "";
                int fileNo = 1;
                String str = "";
                br.readLine();
                br.readLine();
                while((line=br.readLine())!=null)
                {

                    if (line.startsWith("-"))
                    {
                        String outfilename = (String)bwList.get((fileNo-1));
                        FileWriter fileWritter = new FileWriter(outfilename,true);
                        BufferedWriter bw = new BufferedWriter(fileWritter);

                        bw.write(str);
                        bw.close();
                        fileNo = fileNo +1;
                        str = "";
                    }
                    else
                    {
                        str = str + line +"\n";
                    }
                }
                br.close();
                String outfilename = (String)bwList.get((fileNo-1));
                FileWriter fileWritter = new FileWriter(outfilename,true);
                BufferedWriter bw = new BufferedWriter(fileWritter);
                bw.write(str);
                bw.close();
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }
    public HashMap loadFDRMap(String filename)
    {
        HashMap map = new HashMap();

        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pw = splitLine[0];
                double fdr = Double.parseDouble(splitLine[4]);
                map.put(pw, fdr);
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;

    }

    public void writeSigPathwaysFDR_or(double fdrCutoff)
    {

        try
        {
            String filename = DATA_DIR + "bcac_onco_fdr/0.05/bcac-onco-fdr-0.05.txt";
            HashMap map1 = loadFDRMap(filename);
            filename = DATA_DIR + "bcac_icogs_fdr/0.05/bcac-icogs-fdr-0.05.txt";
            HashMap map2 = loadFDRMap(filename);

            Set keys = map1.keySet();
            List pwList1 = new ArrayList(keys);
            String outfilename = DATA_DIR + "bcac_onco_icogs_fdr/0.05/bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-or.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            for (int i = 0; i < pwList1.size();i++)
            {
                String pw = (String)pwList1.get(i);
                if (map2.get(pw)!=null)
                {
                    double fdr1 = (Double)map1.get(pw);
                    double fdr2 = (Double)map2.get(pw);
                    if (fdr1 < fdrCutoff || fdr2 < fdrCutoff)
                    {
                        System.out.println(pw + "\t"+ fdr1 + "\t"+ fdr2);
                        bw.write(pw + "\t"+ fdr1 + "\t"+ fdr2 +"\n");
                    }
                }

            }
            bw.close();

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }
    public List getSigPathwaysFDR_or(double fdrCutoff)
    {
        List pwList = new ArrayList();
        try
        {
            String filename = DATA_DIR + "bcac_onco_icogs_fdr/0.05/bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-or.txt";
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                pwList.add(splitLine[0]);
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return pwList;
    }
    public List getSigPathwaysFDR_and(double fdrCutoff)
    {
        List pwList = new ArrayList();
        try
        {
            String filename = DATA_DIR + "bcac_onco_icogs_fdr/0.05/bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-and.txt";
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                pwList.add(splitLine[0]);
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return pwList;
    }
    public HashMap getSigPathwaystoFDRMap_or(double fdrCutoff)
    {
        HashMap map = new HashMap();
        try
        {
            String filename = DATA_DIR + "bcac_onco_icogs_fdr/0.05/bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-or.txt";

            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pw = splitLine[0];
                double fdr1 = Double.parseDouble(splitLine[1]);
                double fdr2 = Double.parseDouble(splitLine[2]);
                map.put(pw,Math.min(fdr1, fdr2));

            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        Set keys = map.keySet();
        List keyList = new ArrayList(keys);
        System.out.println("Num pathways: " + keyList.size());
        return map;
    }

    public HashMap getSigPathwaystoFDRMap_and(double fdrCutoff)
    {
        HashMap map = new HashMap();
        try
        {
            String filename = DATA_DIR + "bcac_onco_icogs_fdr/0.05/bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-and.txt";

            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pw = splitLine[0];
                double fdr1 = Double.parseDouble(splitLine[1]);
                double fdr2 = Double.parseDouble(splitLine[2]);
                map.put(pw,Math.min(fdr1, fdr2));

            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        Set keys = map.keySet();
        List keyList = new ArrayList(keys);
        System.out.println("Num pathways: " + keyList.size());
        return map;
    }
    public HashMap getMinFDRPerThemeMap(String clusterFile)
    {
        HashMap map = new HashMap();
        try
        {
            List sigPathways = getSigPathwaysFDR_or(0.05);

            HashMap pwToFDRMap = getSigPathwaystoFDRMap_or(0.05);
            System.out.println("Num sig pathways at cutoff: 0.05, " +sigPathways.size());
            HashMap pwToThemeMap = getPathwayToThemeMap(clusterFile);

            Set keys = pwToThemeMap.keySet();
            List keyList = new ArrayList(keys);
            HashMap themeToMinFDRMap = new HashMap();
            for (int i = 0; i < keyList.size();i++)
            {
                String pw = (String)keyList.get(i);
                System.out.println(pw);
                String theme = (String)pwToThemeMap.get(pw);
                if (pw == null)
                    System.out.println("pw is null");
                if (pwToFDRMap ==  null)
                    System.out.println("pwToFDRMap is null");

                Double fdr = (Double)pwToFDRMap.get(pw);
                if (fdr == null)
                    continue;

                Double minFDR = (Double)themeToMinFDRMap.get(theme);
                if (minFDR == null)
                    themeToMinFDRMap.put(theme,fdr);
                else if (fdr < minFDR)
                {
                    themeToMinFDRMap.put(theme,fdr);
                }



            }
            keys = themeToMinFDRMap.keySet();
            keyList = new ArrayList(keys);
            for (int i = 0 ; i < keyList.size();i++)
            {
                String theme = (String)keyList.get(i);
                double minFDR = (Double)themeToMinFDRMap.get(theme);
                //System.out.println(theme + "\t" + minFDR);
                map.put(theme,minFDR);
            }

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;
    }
    public HashMap getPathwayToThemeMap(String clusterThemeFile)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + clusterThemeFile)));
            String line = "";

            String theme = "";
            int numThemes = 0;
            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                {
                    theme = line.substring(1,line.length());
                    numThemes = numThemes+1;
                }
                else
                {
                    if (map.get(line)==null)
                        map.put(line,theme);
                    else
                        System.out.println(line+ " already present in map!");
                }
            }
            br.close();
            System.out.println("Num themes: " + numThemes);

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;
    }

    public void filterPesByFDR(String pesFile, double fdrCutoff)
    {
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + "pes/"+pesFile)));
            String line = "";

            List sigPathways = getSigPathwaysFDR_or(fdrCutoff);

            String outfilename = DATA_DIR + "pes/"+pesFile+".fdr-"+fdrCutoff+"-or.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                double es = Double.parseDouble(splitLine[1]);
                if (es < 0)
                    continue;
                int size = Integer.parseInt(splitLine[2]);
                if (size < 10 || size > 200)
                    continue;
                if (sigPathways.contains(pathway))
                {
                    bw.write(line+"\n");
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
    public HashMap loadPesMap_or(String pesfilename)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR +"pes/"+pesfilename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                double pes = Double.parseDouble(splitLine[1]);
                double size = Integer.parseInt(splitLine[2]);
                map.put(pathway, new double[]{pes, size});
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;
    }
    public HashMap loadFDRMap_or(String fdrORFilename)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + "bcac_onco_icogs_fdr/0.05/"+fdrORFilename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                double fdrOnco = Double.parseDouble(splitLine[1]);
                double fdrIcogs = Double.parseDouble(splitLine[2]);
                map.put(pathway, new double[]{fdrOnco, fdrIcogs});
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;
    }
    public HashMap loadFDRMap(String fdrFilename, String type)
    {
        HashMap map = new HashMap();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + "bcac_"+type+"_fdr/0.05/"+fdrFilename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                double fdr = Double.parseDouble(splitLine[4]);
                map.put(pathway, fdr);
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        return map;
    }
    public void writePesFDRFile_or(String pesFilename, String clusterfile, double fdrCutoff)
    {
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + clusterfile)));
            String line = "";
            String outfilename2 = "bcac-onco-icogs-fdr-0.05-"+fdrCutoff+"-or.txt";
            HashMap map2 = loadFDRMap_or(outfilename2);
            HashMap map3 = loadPesMap_or(pesFilename);
            List sigPathways = getSigPathwaysFDR_or(fdrCutoff);

            String outfilename = DATA_DIR + "overall-pes-fdr-0.05-or.txt";

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            while((line=br.readLine())!=null)
            {
                if (line.startsWith("#"))
                    continue;

                String pathway = line;
                double[] fdrs = (double[]) map2.get(pathway);


                double[] pessize = (double[])map3.get(pathway);
                //System.out.println(pathway);
                double es = pessize[0];
                int size = (int)pessize[1];
                if (sigPathways.contains(pathway))
                {
                    double fdrOnco = fdrs[0];
                    double fdrIcogs = fdrs[1];
                    double minFDR = Math.min(fdrOnco, fdrIcogs);
                    bw.write(pathway + "\t" +es + "\t" + size+"\t"+minFDR+"\n");
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
    public void tableS22(String clusterFile, String pathwayOutputFile, String pesFile)
    {
        try
        {
            HashMap finalMap = new HashMap();

            HashMap map = getSigPathwaystoFDRMap_or(0.05);

            HashMap clusterToPathways = new HashMap();
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + clusterFile)));
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
            br = new BufferedReader(new FileReader(new File(DATA_DIR+ "pes/"+pesFile)));
            line = "";

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathwayName = splitLine[0];
                String pes = splitLine[1];
                pathwayToPESSizeMap.put(pathwayName, pes);
                //System.out.println(pathwayName + ">"+pes);
            }
            br.close();
            br = new BufferedReader(new FileReader(new File(DATA_DIR + pathwayOutputFile)));
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
                    //System.out.println(pathwayName +">*"+line.substring(1,line.length()));
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
                //System.out.println(clusterName0 +">"+pathways0);
                for(int j = 0;j < pathways0.size();j++)
                {
                    String pathway0 = (String)pathways0.get(j);
                    String pathway00 = (String)pathwayToFullName.get(pathway0);
                    if (pathway00 == null)
                    {
                        System.out.println("    pathways null");
                        continue;
                    }
                    String pes00 = (String)pathwayToPESSizeMap.get(pathway0);
                    //System.out.println("*"+pathway0 + ">" + pes00);
                    if (pes00 == null)
                    {
                        //System.out.println("    pes null");
                        continue;
                    }
                    String[] splitpathway00 = pathway00.split("\t");
                    //System.out.println("#Genes: " +splitpathway00[1]);
                    //System.out.println("Score (ES): " +  pes00);
                    //System.out.println("Gene\t#SNPS\tMost Sig. SNP\tP-value\tDistance to Gene (bp)");
                    double fdr = (Double)map.get(pathway0);
                    List genes0 = (List)pathwayToGenes.get(pathway0);
                    //System.out.println("    "+genes0);
                    if (genes0.size()==1)
                    {

                        System.out.println(clusterName0 + "\t" + pathway0 + "\t" + genes0);

                    }
                    else
                    {
                        for (int k = 0 ; k < genes0.size();k++)
                        {
                            //System.out.println(genes0.get(k));
                            String pw = splitpathway00[0];
                            String[] splitPw = pw.split("%");
                            //System.out.println(clusterName0 + "\t" +splitPw[0] + "\t" + splitPw[1] + "\t" + splitPw[2] + "\t" + pes00 + "\t"+ ses + "\t"+ genes0.get(k));
                            List info = (List)finalMap.get(clusterName0);
                            if (info == null)
                                info = new ArrayList();
                            info.add(clusterName0 + "\t" +splitPw[0] + "\t" + splitPw[1] + "\t" + splitPw[2] + "\t" + pes00 + "\t"+ fdr + "\t"+ genes0.get(k));

                            finalMap.put(clusterName0,info);
                        }
                    }
                }

            }
            System.out.println("######");
            keys = finalMap.keySet();
            keyList = new ArrayList(keys);
            Collections.sort(keyList);
            String outfilename = DATA_DIR + "tables/SuppTable22-fdr-0.05-or.txt";
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

    public void tableS23(String clusterFile, String pathwayOutputFile, String pesFile)
    {
        try
        {
            HashMap themeToMinFDRMap = getMinFDRPerThemeMap(clusterFile);
            HashMap clusterToPathways = new HashMap();
            BufferedReader br = new BufferedReader(new FileReader(new File(DATA_DIR + clusterFile)));
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
            br = new BufferedReader(new FileReader(new File(DATA_DIR + "pes/"+pesFile)));
            line = "";

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathwayName = splitLine[0];
                String pes = splitLine[1];
                pathwayToPESSizeMap.put(pathwayName, pes);
            }

            br.close();
            br = new BufferedReader(new FileReader(new File(DATA_DIR + pathwayOutputFile)));
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

                for (int j = 0; j < keyList0.size();j++)
                {
                    double pvalue = (Double)keyList0.get(j);
                    List geneList = (List)valuesToGeneList.get(pvalue);
                    for(int k = 0; k < geneList.size();k++)
                    {
                        List geneList0 = (List)map.get(clusterName0);
                        if (geneList0==null)
                            geneList0 = new ArrayList();
                        geneList0.add(geneList.get(k));
                        map.put(clusterName0,geneList0);

                    }
                }

            }
            String outfilename = DATA_DIR + "tables/SuppTable21-fdr-0.05-or.txt";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));

            keys = map.keySet();
            keyList = new ArrayList(keys);
            Collections.sort(keyList);
            for (int i = 0; i < keyList.size();i++)
            {
                String cluster = (String)keyList.get(i);
                List geneList = (List)map.get(cluster);
                Double fdr = (Double)themeToMinFDRMap.get(cluster);
                if (fdr == null)
                    continue;
                for(int k = 0; k < geneList.size();k++)
                {
                    System.out.println(cluster +"\t"+fdr + "\t"+geneList.get(k));
                    bw.write(cluster +"\t"+fdr + "\t"+geneList.get(k) + "\n");
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
    public void makeEnrichmentFileFDRFilterSingleSigGene(double cutoff)
    {
        try
        {
            HashMap pwToFDRMap = getSigPathwaystoFDRMap_or(cutoff);
            String pathwayOutputFile = DATA_DIR + "overall-predtargets-pathways-siggenes-qcexcluded.txt";
            BufferedReader br = new BufferedReader(new FileReader(new File(pathwayOutputFile)));
            String line = "";
            String pathwayName = "";
            List genes = new ArrayList();
            HashMap pathwayToGenes = new HashMap();
            HashMap pathwayToFullName = new HashMap();
            while((line=br.readLine())!=null)
            {
               // System.out.println(pathwayName +">"+genes);
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
            if (!genes.isEmpty())
            {
                pathwayToGenes.put(pathwayName, genes);
            }
            br.close();

            String filename = DATA_DIR + "pes/all.meta.gwas.icogs.oncoarray.overall.nogenomiccontrol1-HISTONES-5e-2.predicted.targets.txt.pEs_Jan192016_gmt-log2p-qcexcluded.txt.fdr-0.05-or.txt";

            br = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            String outfilename = DATA_DIR + "em/all.meta.gwas.icogs.oncoarray.overall.nogenomiccontrol1-HISTONES-5e-2.predicted.targets.txt.pEs_Jan192016_gmt-log2p-qcexcluded.txt.fdr-"+cutoff+"-or.em-morethanonegene.txt";

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            bw.write("Pathway Id\tDescription\tpvalue\t\t\n");
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String pathway = splitLine[0];
                //System.out.println(pathway);
                List geneList = (List)pathwayToGenes.get(pathway);
                if (geneList.size()==1)
                    continue;
                Double fdr = (Double)pwToFDRMap.get(pathway);
                if (fdr == null)
                    continue;
                bw.write(pathway + "\tNA\t" + fdr+"\n");

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
}
