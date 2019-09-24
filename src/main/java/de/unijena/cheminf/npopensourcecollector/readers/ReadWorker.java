package de.unijena.cheminf.npopensourcecollector.readers;

import java.io.File;
import java.util.ArrayList;

import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.mongocollections.NPDatabase;
import de.unijena.cheminf.npopensourcecollector.mongocollections.NPDatabaseRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.springframework.beans.factory.annotation.Autowired;


public class ReadWorker {


    private File fileToRead;
    private boolean acceptFileFormat = false;
    private String submittedFileFormat ;

    private String fileSource;


    private ArrayList<IAtomContainer> molecules ;

    private Reader reader = null ;


    NPDatabaseRepository npDatabaseRepository;


    public ReadWorker(String fileName){

        npDatabaseRepository = BeanUtil.getBean(NPDatabaseRepository.class);

        this.fileToRead = new File(fileName);

        //System.out.println("\n\n Working on: "+fileToRead.getName() + "\n\n");
        System.out.println("\n\n Working on: "+fileToRead.getAbsolutePath() + "\n\n");

        NPDatabase newDB = new NPDatabase();
        newDB.setLocalFileName(fileToRead.getAbsolutePath());

        npDatabaseRepository.save(newDB);

        acceptFileFormat = acceptFile(fileName);


    }

    public ReadWorker(File file){
        this.fileToRead = file;
        System.out.println("\n\n Working on: "+fileToRead.getAbsolutePath() + "\n\n");

    }




    public boolean startWorker(){
        if(acceptFileFormat){
            return true;
        }
        else{
            return false;
        }

    }

    private boolean acceptFile(String filename) {
        filename = filename.toLowerCase();
        if (filename.endsWith("sdf") || filename.toLowerCase().contains("sdf".toLowerCase())) {
            this.submittedFileFormat="sdf";
            return true;
        } else if (filename.endsWith("smi")  ||
                filename.toLowerCase().contains("smi".toLowerCase()) ||
                filename.toLowerCase().contains("smiles".toLowerCase()) ||
                filename.toLowerCase().contains("smile".toLowerCase())) {
            this.submittedFileFormat="smi";
            return true;
        } else if (filename.endsWith("json")) {
            return false;
        }
        else if (filename.endsWith("mol")  ||
                filename.toLowerCase().contains("mol".toLowerCase())
                || filename.toLowerCase().contains("molfile".toLowerCase())) {
            this.submittedFileFormat="mol";
            return true;
        }
        else if (filename.endsWith("inchi") ||
                filename.toLowerCase().contains("inchi".toLowerCase())
        ){
            this.submittedFileFormat="inchi";
            return true;
        }else if(filename.endsWith("csv") ){
            this.submittedFileFormat="csv";
            return true;
        }


        return false;
    }





    public void doWork(){


        if(this.submittedFileFormat.equals("mol")){
            reader = new MOLReader();
        }
        else if(this.submittedFileFormat.equals("sdf")){
            this.reader = new SDFReader();
        }
        else if(this.submittedFileFormat.equals("smi")){
            reader = new SMILESReader();
        }
        else if(this.submittedFileFormat.equals("inchi")){
            reader = new InChiReader();
        }
        else if(this.submittedFileFormat.equals("csv")){
            reader = new CSVReader();
        }

        this.reader.readFile(this.fileToRead);

    }


    public void doWorkSM(){

        SMReader smReader = new SMReader();
        smReader.readFile(this.fileToRead);

    }


    public String returnSource(){
        return this.reader.returnSource();
    }



}
