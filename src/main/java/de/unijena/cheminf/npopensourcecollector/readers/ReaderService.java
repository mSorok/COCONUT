package de.unijena.cheminf.npopensourcecollector.readers;

import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
public class ReaderService {

    private List<String> molecularFiles;


    public boolean directoryContainsMolecularFiles(String directory){
        boolean molecularFileFound = false;

        try (Stream<Path> walk = Files.walk(Paths.get(directory))) {

            this.molecularFiles = walk.filter(Files::isRegularFile)
                    .map(x -> x.toString()).collect(Collectors.toList());


            for(String f : this.molecularFiles){
                System.out.println(f);
                if(f.contains("sdf") || f.contains("smi") || f.contains("mol") || f.contains("inchi") || f.contains("csv") || f.contains("tsv")){
                    molecularFileFound = true;
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return molecularFileFound;

    }

    public HashSet<String> readMolecularFilesAndInsertInMongo(){

        HashSet<String> totalDatabases = new HashSet<String>();

        ExecutorService taskExecutor = Executors.newFixedThreadPool(this.molecularFiles.size());

        for(String file : this.molecularFiles){
            ReadWorker rw = new ReadWorker();
            rw.setFileToRead(file);

            rw.acceptFileFormat = rw.acceptFile(file);
            boolean start = rw.startWorker();

            if(start){
                taskExecutor.execute(rw);
                System.out.println("Task "+file+" executing");
                String source = file;
                totalDatabases.add(source);
            }

        }

        taskExecutor.shutdown();
        try {
            taskExecutor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        return totalDatabases;

    }

    /*public void readSyntheticMoleculesAndInsertInMongo(String filename){
        //check file extension and if it is not empty
        File smFile = new File(filename);
        if(filename.contains("tsv") && smFile.length()>0){
            //start read worker
            ReadWorker rw = new ReadWorker(smFile);
            rw.doWorkSM();

        }
    }*/


}
