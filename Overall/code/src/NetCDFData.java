import java.io.FileWriter;

import ucar.nc2.*;
import ucar.ma2.Index;
import ucar.ma2.ArrayFloat;
import ucar.ma2.DataType;
import ucar.ma2.ArrayDouble;

import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.BufferedWriter;

/**
 * Reads data stored in NetCDF format.
 * Author: Shirley Hui
 * Date: Sep 20, 2016
 * Time: 1:19:40 PM
 */
public class NetCDFData {

    private final static String DATA_DIR = "/home/shirley/BreastCancer/GSEA/bcancer/";
    private final static String ROOT_DIR = "/home/shirley/BreastCancer/GSEA/Overall/";

    public NetCDFData()
    {

    }
    public List[] getSnpIndexList(String chr, String type)
    {
        List ixList = new ArrayList();
        List snpList = new ArrayList();
        try
        {
            String filename = ROOT_DIR + "data/snps_to_genes_all.meta.gwas.icogs.oncoarray.overall.nogenomiccontrol1-HISTONES-5e-2.predicted.targets.qcsnpsexcluded.txt";
            System.out.println("SNP mapping file: " + filename);

            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            List snpList0 = new ArrayList();

            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String snpId = splitLine[0];
                if (snpId.startsWith(chr+"_"))
                {
                    if (!snpList0.contains(splitLine[0]))
                        snpList0.add(splitLine[0]);
                }
            }
            br.close();
            System.out.println("Num snp ids: " + snpList0.size());
            if (type.equals("onco"))
                filename = DATA_DIR + "bcac_onco/onco_bcac_euro_netcdf_index_chr"+chr+".txt";
            else if (type.equals("icogs"))
                filename = DATA_DIR + "bcac_icogs/icogs_bcac_info_chr"+chr+"_varid.txt";

            br = new BufferedReader(new FileReader(new File(filename)));
            line = "";
            int num =0;
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split(" ");
                String snpId = splitLine[2];
                if (snpList0.contains(snpId))
                {
                    //System.out.println(snpId + "->" + (Integer.parseInt(splitLine[0])-1));
                    ixList.add(Integer.parseInt(splitLine[0])-1);
                    num = num +1;
                    snpList.add(snpId);
                }
            }
            br.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
        System.out.println("Num snp ix: " + ixList.size());
        System.out.println("Num snp id: " + snpList.size());
        return new List[]{ixList, snpList};
    }

    public HashMap getSampleIndexMap(String type)
    {
        HashMap map = new HashMap();
        try
        {
            String filename = ROOT_DIR + "data/sample_orders/bcac_"+type+"_sample_ix_status.txt";
            System.out.println("Sample index file: " + filename);
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line = "";
            while((line=br.readLine())!=null)
            {
                String[] splitLine = line.split("\t");
                String sampleId = splitLine[0];
                int ix = Integer.parseInt(splitLine[1]);
                map.put(ix,sampleId);
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
    public void writeChunk(int incr, String chr, String type, String range, int startIx)
    {
        try
        {
            HashMap sampleMap = getSampleIndexMap(type);
            List[] retList = getSnpIndexList(chr, type);//, range);

            List ixList = retList[0];
            List snpList = retList[1];
            System.out.println("Num samples: " + sampleMap.values().size());
            System.out.println("Num snps: " + snpList.size());
            String filename = DATA_DIR + "bcac_"+type+"/"+type+"_bcac_euro_dosages_chr"+chr+".nc";
            System.out.println("Reading from netCDF file: " + filename);
            System.out.println("Output dir: " + "/mnt/storage/data/bcac_"+type + "/"+range+"-onechunk/");
            NetcdfFile dataFile = NetcdfFile.open(filename, null);
            Variable dataVar = dataFile.findVariable("Dosage");
            if (dataVar == null) {
                System.out.println("Cant find variable Dosage");
                return;
            }
            // Read all the values from the "data" variable into memory.
            int[] shape0 = dataVar.getShape();
            System.out.println(shape0[0] + ", "+ shape0[1]);

            int endIx = startIx +incr;
            int numSnp = 0;
            if (endIx >= ixList.size())
            {
                endIx = ixList.size();
            }
            System.out.println("Chunk: " + startIx + "-" +endIx);

            List ixSubList = new ArrayList(ixList.subList(startIx, endIx));
            List snpSubList = new ArrayList(snpList.subList(startIx, endIx));
            System.out.println("ix sublist size: " + ixSubList.size());
            System.out.println("snp sublist size: " + snpSubList.size());

            int[] shape = new int[2];
            shape[0] = shape0[0];
            shape[1] = 1;

            HashMap map = new HashMap();
            int num = 0;

            for (int i = 0; i < ixSubList.size();i++)
            {
                int ix = (Integer)ixSubList.get(i);

                int[] origin = new int[2];
                origin[0] = 0;
                origin[1] = ix;
                ArrayFloat.D2 dataArray = (ArrayFloat.D2) dataVar.read(origin, shape);
                num= num +1;
                int numSamples =0;

                for (int j = 0; j < shape0[0]; j++)
                {
                    int ixx = j+1;
                    String sampleId = (String)sampleMap.get(ixx);
                    if (sampleId == null)
                        continue;

                    numSamples = numSamples+1;
                    double value = dataArray.get(j,0);
                    List list = (List)map.get(ixx);
                    if (list==null)
                        list = new ArrayList();
                    list.add(value);
                    map.put(ixx,list);
                }
                numSnp = numSnp+1;
                System.out.println(numSnp + ": " +i +": "+ snpSubList.get(i) + ": origin: ["+origin[0] + ","+origin[1]+"], shape: [" + shape[0] + "," + shape[1] +"], num samples: " + numSamples);

            }
            System.out.println("Num snps: " + num);
            String outfilename = "/mnt/storage/data/bcac_"+type + "/"+range+"-onechunk/"+type+"_bcac_euro_dosages_chr"+chr+"_"+range+"-"+endIx+".txt";

            System.out.println("Writing to file: " + outfilename);

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            System.out.println("Num header tokens: " +(snpSubList.size()+1));
            bw.write("PID");
            for (int i= 0 ; i < snpSubList.size();i++)
            {
                bw.write("\t"+snpSubList.get(i));
            }
            bw.write("\n");
            Set keys = map.keySet();
            List keyList = new ArrayList(keys);
            Collections.sort(keyList);
            for (int i=0;i < keyList.size();i++)
            {
                int ix = (Integer)keyList.get(i);
                String sampleId = (String)sampleMap.get(ix);
                if (sampleId == null)
                    continue;

                List list = (List)map.get(ix);
                //System.out.println("ix="+ix + ", sampleid="+sampleId+", size=" +list.size());
                bw.write(sampleId);
                for (int j=0;j < list.size();j++)
                {
                    bw.write("\t"+list.get(j));
                }
                bw.write("\n");
            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public void writeAllChunks(int incr, String chr, String type, String range)
    {
        try
        {
            HashMap sampleMap = getSampleIndexMap(type);
            List[] retList = getSnpIndexList(chr, type);//, range);

            List ixList = retList[0];
            List snpList = retList[1];
            System.out.println("Num samples: " + sampleMap.values().size());
            System.out.println("Num snps: " + snpList.size());

            String filename = DATA_DIR + "bcac_"+type+"/"+type+"_bcac_euro_dosages_chr"+chr+".nc";
            System.out.println("Reading from netCDF file: " + filename);
            NetcdfFile dataFile = NetcdfFile.open(filename, null);
            Variable dataVar = dataFile.findVariable("Dosage");
            if (dataVar == null) {
                System.out.println("Cant find variable Dosage");
                return;
            }
            // Read all the values from the "data" variable into memory.
            int[] shape0 = dataVar.getShape();
            System.out.println(shape0[0] + ", "+ shape0[1]);

            int startIx = 0;
            int endIx = startIx +incr;
            boolean stop = false;
            int numSnp = 0;
            while(!stop)
            {
                if (endIx >= ixList.size())
                {

                    endIx = ixList.size();
                    stop = true;
                }
                System.out.println("Chunk: " + startIx + "-" +endIx);

                List ixSubList = new ArrayList(ixList.subList(startIx, endIx));
                List snpSubList = new ArrayList(snpList.subList(startIx, endIx));
                System.out.println("ix sublist size: " + ixSubList.size());
                System.out.println("snp sublist size: " + snpSubList.size());

                int[] shape = new int[2];
                shape[0] = shape0[0];
                shape[1] = 1;

                HashMap map = new HashMap();
                int num = 0;
                for (int i = 0; i < ixSubList.size();i++)
                {
                    int ix = (Integer)ixSubList.get(i);

                    int[] origin = new int[2];
                    origin[0] = 0;
                    origin[1] = ix;
                    ArrayFloat.D2 dataArray = (ArrayFloat.D2) dataVar.read(origin, shape);
                    num= num +1;
                    int numSamples =0;
                    for (int j = 0; j < shape0[0]; j++)
                    {
                        int ixx = j+1;
                        String sampleId = (String)sampleMap.get(ixx);
                        if (sampleId == null)
                            continue;


                        numSamples = numSamples+1;
                        double value = dataArray.get(j,0);
                        List list = (List)map.get(ixx);
                        if (list==null)
                            list = new ArrayList();
                        list.add(value);
                        map.put(ixx,list);
                    }
                    numSnp = numSnp+1;
                    System.out.println(numSnp + ": " +i +": "+ snpSubList.get(i) + ": origin: ["+origin[0] + ","+origin[1]+"], shape: [" + shape[0] + "," + shape[1] +"], num samples: " + numSamples);

                }
                System.out.println("Num snps: " + num);

                String outfilename = "/mnt/storage/data/bcac_"+type + "/"+range+"/"+type+"_bcac_euro_dosages_chr"+chr+"_"+range+"-"+endIx+".txt";

                System.out.println("Writing to file: " + outfilename);

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
                System.out.println("Num header tokens: " +(snpSubList.size()+1));
                bw.write("PID");
                for (int i= 0 ; i < snpSubList.size();i++)
                {
                    bw.write("\t"+snpSubList.get(i));
                }
                bw.write("\n");
                Set keys = map.keySet();
                List keyList = new ArrayList(keys);
                Collections.sort(keyList);
                for (int i=0;i < keyList.size();i++)
                {
                    int ix = (Integer)keyList.get(i);
                    String sampleId = (String)sampleMap.get(ix);
                    if (sampleId == null)
                        continue;

                    List list = (List)map.get(ix);
                    //System.out.println("ix="+ix + ", sampleid="+sampleId+", size=" +list.size());
                    bw.write(sampleId);
                    for (int j=0;j < list.size();j++)
                    {
                        bw.write("\t"+list.get(j));
                    }
                    bw.write("\n");
                }
                bw.close();
                startIx = endIx;
                endIx = startIx + incr;

            }


        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public void readTest()
    {
        try
        {
            String filename = DATA_DIR + "bcac_icogs/icogs_bcac_euro_dosages_chr1.nc";
            NetcdfFile dataFile = NetcdfFile.open(filename, null);
            Variable dataVar = dataFile.findVariable("Dosage");
            if (dataVar == null) {
                System.out.println("Cant find variable Dosage");
                return;
            }
            // Read all the values from the "data" variable into memory.
            int[] shape0 = dataVar.getShape();
            System.out.println(shape0[0] + "," + shape0[1]);
            System.out.println(dataVar.toString()) ;

        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

    public static void main(String[] args)
    {
        NetCDFData n = new NetCDFData();
        int incr = Integer.parseInt(args[0]);
        String chr = args[1];
        String type = args[2];//icogs or onco
        String range = "0.05";//args[3];//0.001,0.001-0.005
        int startIx = Integer.parseInt(args[3]);
        System.out.println("Writing chunks of size: " + incr);
        n.writeChunk(incr,chr,type, range, startIx);
        //n.writeAllChunks(incr,chr, type,range);
        //n.readTest();
    }
}

