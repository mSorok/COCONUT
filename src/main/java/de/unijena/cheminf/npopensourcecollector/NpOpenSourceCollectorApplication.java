package de.unijena.cheminf.npopensourcecollector;

import com.mongodb.MongoClientOptions;
import de.unijena.cheminf.npopensourcecollector.readers.ReaderService;
import de.unijena.cheminf.npopensourcecollector.services.*;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;


import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
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


    @Autowired
    ExportService exportService;

    @Autowired
    AnnotationLevelService annotationLevelService;

    @Autowired
    NamingService namingService;


    @Autowired
    FingerprintsCountsFiller fingerprintsCountsFiller;


    @Autowired
    IntegratePassPredictions integratePassPredictions;



    public static void main(String[] args) {
        SpringApplication.run(NpOpenSourceCollectorApplication.class, args);
    }

    @Override
    public void run(String... args) throws Exception {





        System.out.println("Code version from 15th of March 2021");

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
            else if(args[0].equals("fetchPASSpredictions")){
                integratePassPredictions.readAndInsertPassPredictions(args[1]);
            }
            else if(args[0].equals("cleanRecomputeMissing")){


                System.out.println("Fragmenting everything from scratch");
                fragmentCalculatorService.doParallelizedWork(42);

                System.out.println("Done fragmenting natural products");
                SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
                System.out.println("at: "+formatter.format(new Date())+"\n");


                //evaluate annotation level
                System.out.println("evaluating annotation levels");
                annotationLevelService.doWorkForAll();
                System.out.println("done");

                // Compute additional features
                molecularFeaturesComputationService.doWork();
                updaterService.updateSourceNaturalProductsParallelized(42);


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
            else if(args[0].equals("updateCNPid")){
                System.out.println("Updating COCONUT ids");
                createCNPidService.clearIDs();
                createCNPidService.importIDs("coconut_ids_june2020.csv");
                createCNPidService.createIDforNewMolecules();

                updaterService.updateSourceNaturalProductsParallelized(40);
                while (!updaterService.processFinished()) {
                    System.out.println("I'm waiting");
                    TimeUnit.MINUTES.sleep(1);
                }
                System.out.println("done");
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
                //readerService.readSyntheticMoleculesAndInsertInMongo(args[1]); //tsv file
                //molecularFeaturesComputationService.doWorkForSM();
                System.out.println("Need to reimplement");
            }
            else if(args[0].equals("generateSDF")){

                exportService.generateSDF("COCONUT.sdf");

            }
            else if(args[0].equals("exportAllFormats")){
                exportService.generateSimpleSMILESFile("COCONUT.smi");
                exportService.generateSDF("COCONUT.sdf");
                exportService.generateFullSMILESfile("COCONUTstereo.smi");
            }
            else if(args[0].equals("generateTSV")){
                exportService.generateTSV("COCONUT.tsv");
            }
            else if(args[0].equals("onlyImportCoconutIds")){
                int index_of_id_file = Arrays.asList(args).indexOf("onlyImportCoconutIds")+1;
                createCNPidService.importIDs(args[index_of_id_file]);
                createCNPidService.createIDforNewMolecules();

            }
            else if(args[0].equals("updateBitFingerprints")){

                molecularFeaturesComputationService.convertToBitSet();

            }else if(args[0].equals("createPubchemBitCounts")) {
                molecularFeaturesComputationService.createPubchemBitCounts();

            }else if (args[0].equals("evaluateAnnotation")){
                annotationLevelService.doWorkForAll();
            }else if (args[0].equals("namesToLowerCase")){
                namingService.namesToLowcase();
            }
            /*else if(args[0].equals("generateUniqueSmiles")){
                molecularFeaturesComputationService.generateUniqueSmiles();
            }*/
            else { //Filling from scratch
                //cleaning the DB before filling it



                String dataDirectory = args[0];

                boolean canContinue = readerService.directoryContainsMolecularFiles(dataDirectory);


                if (canContinue) {
                    //insert in mongodb


                    mongoTemplate.getDb().drop();

                    readerService.readMolecularFilesAndInsertInMongo();

                    //unify
                    //npUnificationService.fetchSourceNames();
                    npUnificationService.doWork();

                    //molecularFeaturesComputationService.generateUniqueSmiles();


                    fragmentReaderService.doWork(0, args[2]);
                    fragmentReaderService.doWork(1, args[3]);


                    //fragmentCalculatorService.doWork();

                    fragmentCalculatorService.doParallelizedWork(40);


                    System.out.println("Done fragmenting natural products");
                    SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
                    System.out.println("at: "+formatter.format(new Date())+"\n");

                    if(Arrays.asList(args).contains("importCOCONUTids")) {
                        //coconut_ids_april2020.csv
                        System.out.println("importing  old COCONUT ids");
                        int index_of_id_file = Arrays.asList(args).indexOf("importCOCONUTids")+1;
                        createCNPidService.clearIDs();
                        createCNPidService.importIDs(args[index_of_id_file]);
                        createCNPidService.createIDforNewMolecules();

                        System.out.println("done importing IDS and generating news ones");
                    }else{

                        if(Arrays.asList(args).contains("idPrefix")){
                            int index_of_prefix = Arrays.asList(args).indexOf("idPrefix")+1;
                            createCNPidService.createDeNovoIDs(args[index_of_prefix]);
                        }else {

                            createCNPidService.createDeNovoIDs();
                        }
                    }

                    //evaluate annotation level
                    System.out.println("evaluating annotation levels");
                    annotationLevelService.doWorkForAll();
                    System.out.println("done");

                    // Compute additional features
                    molecularFeaturesComputationService.doWork();
                    updaterService.updateSourceNaturalProductsParallelized(42);

                    namingService.namesToLowcase();

                    fingerprintsCountsFiller.doWork();





                    molecularFeaturesComputationService.createPubchemBitCounts();

                    if(Arrays.asList(args).contains("passPredictions")){
                        int index_of_passfile = Arrays.asList(args).indexOf("passPredictions")+1;
                        integratePassPredictions.readAndInsertPassPredictions(args[index_of_passfile]);
                    }




                    //read and insert synthetic molecules
                    //readerService.readSyntheticMoleculesAndInsertInMongo(args[1]); //tsv file
                    //molecularFeaturesComputationService.doWorkForSM();

                    //molecularFeaturesComputationService.convertToBitSet();

/*
                    //compute similarities between natural products
                    similarityComputationService.generateAllPairs();
                    // //similarityComputationService.computeSimilarities();
                    similarityComputationService.doParallelizedWork(40);
                    while (!similarityComputationService.processFinished()) {
                        System.out.println("I'm waiting");
                        TimeUnit.MINUTES.sleep(1);
                    }
*/

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
