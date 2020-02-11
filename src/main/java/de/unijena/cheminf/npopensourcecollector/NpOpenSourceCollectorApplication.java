package de.unijena.cheminf.npopensourcecollector;

import com.mongodb.MongoClientOptions;
import de.unijena.cheminf.npopensourcecollector.readers.ReaderService;
import de.unijena.cheminf.npopensourcecollector.services.*;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;


import java.util.concurrent.TimeUnit;

import static java.lang.System.exit;

@SpringBootApplication
public class NpOpenSourceCollectorApplication implements CommandLineRunner {

    @Autowired
    ReaderService readerService;


    @Autowired
    NPUnificationService npUnificationService;

    @Autowired
    FragmentReaderService fragmentReaderService;

    @Autowired
    FragmentCalculatorService fragmentCalculatorService;

    @Autowired
    MolecularFeaturesComputationService molecularFeaturesComputationService;

    @Autowired
    org.springframework.data.mongodb.core.MongoTemplate mongoTemplate;

    @Autowired
    SimilarityComputationService similarityComputationService;

    @Autowired
    UpdaterService updaterService;

    @Autowired
    CreateCNPidService createCNPidService;

    public static void main(String[] args) {
        SpringApplication.run(NpOpenSourceCollectorApplication.class, args);
    }

    @Override
    public void run(String... args) throws Exception {





        System.out.println("Code version from 8th October 2019");

        if (args.length > 0) {

            if(args[0].equals("recomputeMissing")){
                fragmentCalculatorService.doWorkRecompute();
                molecularFeaturesComputationService.doWorkRecompute();
                updaterService.updateSourceNaturalProductsParallelized(40);
                while (!updaterService.processFinished()) {
                    System.out.println("I'm waiting");
                    TimeUnit.MINUTES.sleep(1);
                }

            }
            else if(args[0].equals("addCNPid")){

                System.out.println("Creating de novo COCONUT IDs");
                createCNPidService.createDeNovoIDs();
                updaterService.updateSourceNaturalProductsParallelized(40);
                while (!updaterService.processFinished()) {
                    System.out.println("I'm waiting");
                    TimeUnit.MINUTES.sleep(1);
                }

            }
            else if(args[0].equals("runOnlySimilarity")){
                //compute similarities between natural products
                similarityComputationService.generateAllPairs();
                // //similarityComputationService.computeSimilarities();
                similarityComputationService.doParallelizedWork(40);
                while (!similarityComputationService.processFinished()) {
                    System.out.println("I'm waiting");
                    TimeUnit.MINUTES.sleep(1);
                }
            }
            else if(args[0].equals("onlyAddSM")){
                //read and insert synthetic molecules
                readerService.readSyntheticMoleculesAndInsertInMongo(args[1]); //tsv file
                molecularFeaturesComputationService.doWorkForSM();
            }
            else { //Filling from scratch
                //cleaning the DB before filling it
                mongoTemplate.getDb().drop();

                String dataDirectory = args[0];

                boolean canContinue = readerService.directoryContainsMolecularFiles(dataDirectory);


                if (canContinue) {
                    //insert in mongodb


                    readerService.readMolecularFilesAndInsertInMongo();

                    //unify
                    //npUnificationService.fetchSourceNames();
                    npUnificationService.doWork();


                    fragmentReaderService.doWork(0, args[2]);
                    fragmentReaderService.doWork(1, args[3]);


                    //fragmentCalculatorService.doWork();

                    fragmentCalculatorService.doParallelizedWork(40);

                    while (!fragmentCalculatorService.processFinished()) {
                        System.out.println("I'm waiting");
                        TimeUnit.MINUTES.sleep(1);
                    }


                    molecularFeaturesComputationService.doWork();
                    createCNPidService.createDeNovoIDs();
                    updaterService.updateSourceNaturalProductsParallelized(40);

                    //read and insert synthetic molecules
                    readerService.readSyntheticMoleculesAndInsertInMongo(args[1]); //tsv file
                    molecularFeaturesComputationService.doWorkForSM();


                    //compute similarities between natural products
                    similarityComputationService.generateAllPairs();
                    // //similarityComputationService.computeSimilarities();
                    similarityComputationService.doParallelizedWork(40);
                    while (!similarityComputationService.processFinished()) {
                        System.out.println("I'm waiting");
                        TimeUnit.MINUTES.sleep(1);
                    }

                } else {
                    System.out.println("Could not find files with molecules in the provided directory!");
                    exit(0);
                }

            }


            System.out.println("Normal exit");
            exit(0);

        }




    }
}
