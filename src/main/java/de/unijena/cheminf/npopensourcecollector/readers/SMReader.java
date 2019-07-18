package de.unijena.cheminf.npopensourcecollector.readers;

import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeChecker;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SyntheticMolecule;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SyntheticMoleculeRepository;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.qsar.result.IntegerResult;

import java.io.*;

public class SMReader {

    File file;
    MoleculeChecker moleculeChecker;
    SyntheticMoleculeRepository syntheticMoleculeRepository;

    private LineNumberReader tsvReader;


    public SMReader(){
        moleculeChecker = BeanUtil.getBean(MoleculeChecker.class);
        syntheticMoleculeRepository =  BeanUtil.getBean(SyntheticMoleculeRepository.class);

    }

    public void readFile(File file){

        this.file = file;
        int count = 1;
        String line;


        try{
            tsvReader = new LineNumberReader(new InputStreamReader(new FileInputStream(file)));
            System.out.println("Reading and inserting synthetic molecules from ZINC");


            while ((line = tsvReader.readLine()) != null  && count <= 600000) {

                String[] splitted = line.split("\\s+");
                //smiles, inchi, inchikey, contains_sugar, total_atom_number, heavy_atom_number, sugar_free_total_atom_number, sugar_free_heavy_atom_number, npl_score, npl_sugar_score, npl_noh_score, molecular_weight, molecular_formula, number_of_carbons, number_of_nitrogens, number_of_oxygens, number_of_rings, number_repeated_fragments

                SyntheticMolecule newSM = new SyntheticMolecule();

                newSM.setSmiles(splitted[0]);
                newSM.setInchi(splitted[1]);
                newSM.setInchikey(splitted[2]);
                newSM.setContains_sugar(Integer.parseInt(splitted[3]));
                newSM.setTotal_atom_number(Integer.parseInt(splitted[4]));
                newSM.setHeavy_atom_number(Integer.parseInt(splitted[5]));
                newSM.setSugar_free_total_atom_number(Integer.parseInt(splitted[6]));
                newSM.setSugar_free_heavy_atom_number(Integer.parseInt(splitted[7]));
                newSM.setNpl_score(Double.parseDouble(splitted[8]));
                newSM.setNpl_sugar_score(Double.parseDouble(splitted[9]));
                newSM.setNpl_noh_score(Double.parseDouble(splitted[10]));
                newSM.setMolecular_weight(Double.parseDouble(splitted[11]));
                newSM.setMolecular_formula(splitted[12]);
                newSM.setNumber_of_carbons(Integer.parseInt(splitted[13]));
                newSM.setNumber_of_nitrogens(Integer.parseInt(splitted[14]));
                newSM.setNumber_of_oxygens(Integer.parseInt(splitted[15]));
                newSM.setNumber_of_rings(Integer.parseInt(splitted[16]));
                newSM.setNumber_repeated_fragments(Integer.parseInt(splitted[17]));

                syntheticMoleculeRepository.save(newSM);

            }



        } catch (IOException ex) {
            System.out.println("Oops ! File not found. Please check if the -in file or -out directory is correct");
            ex.printStackTrace();
            System.exit(0);

        }






    }
}
