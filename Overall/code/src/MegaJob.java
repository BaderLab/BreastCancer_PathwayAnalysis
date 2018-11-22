import java.io.*;
import java.util.*;

/**
 * Code to make files used to run batch jobs on Scinet
 * Author: shirleyhui
 * Date: Jul 22, 2016
 * Time: 1:18:58 PM
 */
public class MegaJob {

    private static String DATA_DIR = "/scratch/g/gbader/shui/data/";
    private static String ROOT_DIR = "/scratch/g/gbader/shui/";
    public MegaJob()
    {

    }
    public static void main(String[] args)
    {
        MegaJob j = new MegaJob();
        String type = args[0];
	String range = args[1];
        //j.make_chisq(type,range);
        //j.submit_chisq(type,range); 
	
	// For 1000 perms
        int startIx = 1;
	for (int i = 0; i < 200;i++)
        {
	    j.make_gsea(startIx, type, range);
            startIx = startIx + 5;
        }
		
        j.make_submit_script(1,1, type, range);
        j.make_submit_script(251,2, type, range);
        j.make_submit_script(501,3, type, range);
        j.make_submit_script(751,4, type, range);
         
  }
    public void make_chisq(String type,String range)
    {
        try
        {
            String outfilename = "";
            BufferedWriter bw = null;

            int written = 8;
            String str = "";
            int fileno = 1;

            File dirFile = new File(DATA_DIR + "bcac_"+type+"/"+range+"/");
            File[] files = dirFile.listFiles();

            for (int i = 0;i < files.length;i++)
            {

                if (written == 8)
                {
                    if (!str.equals(""))
                    {
                        bw.write(str);
                        bw.write("EOF");
                        bw.close();
                    }
                    outfilename = ROOT_DIR + "job_scripts/bcac_"+type+"_chisq/"+range+"/run-chisq-"+range+"-"+fileno + ".sh";
                    bw = new BufferedWriter(new FileWriter(new File(outfilename)));
                    bw.write("#!/bin/bash\n");
                    bw.write("cd $PBS_O_WORKDIR\n");
                    bw.write("export OMP_NUM_THREADS=1\n");
                    bw.write("module load intel/15.0.6\n");
		    bw.write("module load curl/7.49.1\n");
                    bw.write("module load gnu-parallel/20140622\n");
                    bw.write("module load R/3.3.0\n");
                    bw.write("module load extras\n");      
                    bw.write("parallel -j 8 <<EOF\n");
                    written = 0;
                    str = "";
                    fileno = fileno+1;
                }

                String filename = files[i].getName();
                //System.out.println("*"+filename);
		String[] splitName = filename.split("_");
                String chr = splitName[4];
		//System.out.println("**"+splitName[5]);
	 	splitName = splitName[5].split("-");
		//System.out.println(splitName[2]);
		System.out.println("***"+splitName[1]);
                String incr = splitName[1].substring(0,splitName[1].length()-4);  
                //String incr = splitName[2].substring(0,splitName[2].length()-4); 
                String suffix = chr+"-"+incr;
                str = str +"mkdir ChiSq-"+suffix+"; cd ChiSq-" +suffix +"; Rscript ../R_ChiSq-1000.R ../data/bcac_"+type+"/"+range+"/"+filename+"; echo \"job "+suffix+" finished\"\n";

                written = written +1;

            }
            bw.write(str);
            bw.write("EOF");
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }
    public void submit_chisq(String type, String range)
    {
	try
	{
   	   String outfilename = ROOT_DIR + "submit-bcac-"+type+"-chisq-"+range+".sh";
	   BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
	   File dirFile = new File(ROOT_DIR + "job_scripts/bcac_"+type+"_chisq/"+range+"/");
           File[] files = dirFile.listFiles();

           for (int i = 0;i < files.length;i++)
           {
	      String file= files[i].getName();
	      bw.write("qsub -l nodes=1:ppn=8,walltime=4:30:00 -o oo-"+i+".txt -e ee-"+i+".txt job_scripts/bcac_"+type+"_chisq/"+range+"/"+ file +"\n");

	   }
	   bw.close();
	}
	catch(Exception e)
	{
	   System.out.println("Exception: " + e);
	   e.printStackTrace();
	}


    }
    public void make_gsea(int startIx, String type, String range)
    {
        try
        {


            String outfilename =ROOT_DIR+ "job_scripts/bcac_"+type+"_gsea/"+range+"/run-gsea-"+range+"-"+startIx + ".sh";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            bw.write("#!/bin/bash\n");
            bw.write("cd $PBS_O_WORKDIR\n");
            bw.write("export OMP_NUM_THREADS=1\n");
            bw.write("module load intel\n");
            bw.write("module load openmpi\n");
            bw.write("module load gcc\n");
            bw.write("module load gnu-parallel/20140622\n");
            bw.write("module load python\n");
            bw.write("parallel -j 1 <<EOF\n");

            bw.write("mkdir Gsea0; cd Gsea0; python ../mygseaCalcESCentrality-perm.py ../bcac_"+type+"_chisq_"+range+"/perm/bcac_"+type+"_euro_dosages_"+range+".perm0.txt; echo \"job 0 finished\"\n");
            for (int i = 0;i<5;i++)
            {
                bw.write("mkdir Gsea"+(startIx+i)+"; cd Gsea" +(startIx+i) +"; python ../mygseaCalcESCentrality-perm.py ../bcac_"+type+"_chisq_"+range+"/perm/bcac_"+type+"_euro_dosages_"+range+".perm"+(startIx+i)+".txt; echo \"job "+(startIx+i)+" finished\"\n");
            }

            bw.write("EOF");
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }
    }

   public void make_submit_script(int ix, int num,String type, String range)
    {
        try
        {

            String outfilename = ROOT_DIR +"submit-bcac-"+type+"-gsea-"+range+"-"+num+".sh";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfilename)));
            for (int i = 0; i < 50;i++)
            {
                bw.write("qsub -l nodes=1:ppn=8,walltime=4:00:00 -o oo-"+ix+".txt -e ee-"+ix+".txt job_scripts/bcac_"+type+"_gsea/"+range+"/run-gsea-"+range+"-"+ix+".sh\n");
                ix  = ix + 5;
            }
            bw.close();
        }
        catch(Exception e)
        {
            System.out.println("Exception: " + e);
            e.printStackTrace();
        }

    }

}
