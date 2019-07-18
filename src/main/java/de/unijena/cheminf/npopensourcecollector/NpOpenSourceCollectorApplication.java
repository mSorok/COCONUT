package de.unijena.cheminf.npopensourcecollector;

import de.unijena.cheminf.npopensourcecollector.readers.ReaderService;
import de.unijena.cheminf.npopensourcecollector.services.FragmentCalculatorService;
import de.unijena.cheminf.npopensourcecollector.services.FragmentReaderService;
import de.unijena.cheminf.npopensourcecollector.services.MolecularFeaturesComputationService;
import de.unijena.cheminf.npopensourcecollector.services.NPUnificationService;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;



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

    public static void main(String[] args) {
        SpringApplication.run(NpOpenSourceCollectorApplication.class, args);
    }

    @Override
    public void run(String... args) throws Exception {

        //cleaning the DB before filling it
       // mongoTemplate.getDb().drop();

        System.out.println("Code version from 21 June 2019");

        if (args.length > 0) {
            String dataDirectory = args[0];

            boolean canContinue = readerService.directoryContainsMolecularFiles(dataDirectory);


            if(canContinue){
                //insert in mongodb
                readerService.readMolecularFilesAndInsertInMongo();

                //unify
                npUnificationService.doWork();


                fragmentReaderService.doWork(0, "/Users/maria/Projects/NPOpenSourceCollector/fragments/fragment_without_sugar.txt");
                fragmentReaderService.doWork(1, "/Users/maria/Projects/NPOpenSourceCollector/fragments/fragment_with_sugar.txt");


                fragmentCalculatorService.doWork();

                molecularFeaturesComputationService.doWork();


                //read and insert synthetic molecules
                readerService.readSyntheticMoleculesAndInsertInMongo(args[1]); //tsv file
                molecularFeaturesComputationService.doWorkForSM();


            }
            else{
                System.out.println("Could not find files with molecules in the provided directory!");
                exit(0);
            }




            System.out.println("Normal exit");
            exit(0);

        }




    }
}
