package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.convert.ConversionFailedException;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


@Service
public class ExportService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;


    public void generateTSV(String filename){
        System.out.println("Generating TSV file for all COCONUT");

        FileWriter fw = null;
        try {
            fw = new FileWriter(filename);

            List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();
            int count = 0;
            for(UniqueNaturalProduct np : allNP){
                fw.write(np.coconut_id+"\t"+np.clean_smiles+"\n");
            }

            fw.close();
            System.out.println("Successfully wrote to the file.");


        } catch (IOException e) {
            e.printStackTrace();
        }




        System.out.println("done");

    }

    public void generateFullSDF(String filename){

        System.out.println("Generating the full SDF file for all stereo COCONUT");

        //TODO

    }

    public void generateFullSMILESfile(String filename){

        System.out.println("Generating all stereo SMILES file for all COCONUT");



        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

            List<String> coconut_ids = uniqueNaturalProductRepository.findAllCoconutIds();
            int count = 0;

            for(String coconut_id : coconut_ids){
                UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);

                String baseSmiles = np.unique_smiles;

                writer.write(baseSmiles+" "+coconut_id+"\n");

                int counter=1;

                for(String absoluteSmiles : np.absolute_smiles.keySet()){
                    if(!np.absolute_smiles.get(absoluteSmiles).equals("nostereo") && !np.absolute_smiles.get(absoluteSmiles).equals(np.unique_smiles)){

                        writer.write(np.absolute_smiles.get(absoluteSmiles) + " "+coconut_id+"."+counter + "\n");
                        counter = counter+1;
                    }


                }

                count++;

                if (count % 50000 == 0) {
                    System.out.println("Molecules written: " + count);
                }

            }

            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void generateSimpleSMILESFile(String filename){

        System.out.println("Generating SMILES file for all COCONUT");



        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

            List<String> coconut_ids = uniqueNaturalProductRepository.findAllCoconutIds();
            int count = 0;

            for(String coconut_id : coconut_ids){
                UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);

                writer.write(np.unique_smiles+" "+coconut_id+"\n");

                count++;

                if (count % 50000 == 0) {
                    System.out.println("Molecules written: " + count);
                }

            }

            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }


    }


    public void generateSDF(String filename){
        System.out.println("Generating SDF file for all COCONUT");



        FileWriter fw = null;
        try {
            fw = new FileWriter(filename);

            SDFWriter writer = new SDFWriter(fw);

            //List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();

            List<String> coconut_ids = uniqueNaturalProductRepository.findAllCoconutIds();
            int count = 0;

            for(String coconut_id : coconut_ids){
                try {
                    UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);
                    IAtomContainer ac = atomContainerToUniqueNaturalProductService.createAtomContainer(np);

                    // add most of molecular descriptors and available metadata
                    ac.setProperty("coconut_id", np.coconut_id);
                    ac.setProperty("inchi", np.inchi);
                    ac.setProperty("inchikey", np.inchikey);
                    ac.setProperty("SMILES", np.clean_smiles);
                    ac.setProperty("sugar_free_smiles", np.sugar_free_smiles);
                    ac.setProperty("molecular_formula", np.molecular_formula);
                    ac.setProperty("molecular_weight", np.molecular_weight);
                    ac.setProperty("citationDOI", np.citationDOI.toString());
                    ac.setProperty("textTaxa", np.textTaxa.toString());
                    ac.setProperty("geoLocation", np.geoLocation);
                    ac.setProperty("name", np.name);
                    ac.setProperty("synonyms", np.synonyms.toString());
                    ac.setProperty("NPL score", np.npl_score);
                    ac.setProperty("number_of_carbons", np.number_of_carbons);
                    ac.setProperty("number_of_nitrogens", np.number_of_nitrogens);
                    ac.setProperty("number_of_oxygens", np.number_of_oxygens);
                    ac.setProperty("number_of_rings", np.number_of_rings);
                    ac.setProperty("total_atom_number", np.total_atom_number);
                    ac.setProperty("bond_count", np.bond_count);
                    ac.setProperty("found_in_databases", np.found_in_databases.toString());
                    ac.setProperty("murko_framework", np.murko_framework);
                    ac.setProperty("alogp", np.alogp);
                    ac.setProperty("apol", np.apol);
                    ac.setProperty("topoPSA", np.topoPSA);

                    writer.write(ac);

                    count++;

                    if (count % 50000 == 0) {
                        System.out.println("Molecules written: " + count);
                    }
                }catch(ConversionFailedException e){
                    System.out.println(coconut_id);
                }


            }

            writer.close();

        } catch (IOException | CDKException e) {
            e.printStackTrace();
        }




    System.out.println("done");

    }
}
