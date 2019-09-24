package de.unijena.cheminf.npopensourcecollector.readers;

import de.unijena.cheminf.npopensourcecollector.misc.DatabaseTypeChecker;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.services.AtomContainerToSourceNaturalProductService;
import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeChecker;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import net.sf.jniinchi.INCHI_OPTION;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.springframework.context.annotation.Bean;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.*;

public class SDFReader implements Reader{

    File file;
    ArrayList<IAtomContainer> listOfMolecules;

    private IteratingSDFReader reader = null;

    SourceNaturalProductRepository sourceNaturalProductRepository;
    AtomContainerToSourceNaturalProductService ac2snp;


    MoleculeChecker moleculeChecker;

    DatabaseTypeChecker databaseTypeChecker;

    String source;




    public SDFReader(){

        this.listOfMolecules = new ArrayList<IAtomContainer>();
        sourceNaturalProductRepository = BeanUtil.getBean(SourceNaturalProductRepository.class);
        ac2snp = BeanUtil.getBean(AtomContainerToSourceNaturalProductService.class);
        moleculeChecker = BeanUtil.getBean(MoleculeChecker.class);
        databaseTypeChecker = BeanUtil.getBean(DatabaseTypeChecker.class);


    }

    @Override
    public void readFile(File file) {

        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique); //Unique - canonical SMILES string, different atom ordering produces the same* SMILES. No isotope or stereochemistry encoded.

        this.file = file;
        int count = 1;

        this.source = file.getName().toLowerCase().replace(".sdf", "");


        try{

            reader = new IteratingSDFReader(new FileInputStream(file), DefaultChemObjectBuilder.getInstance());
            System.out.println("SDF reader creation and inserting in MongoDB for "+source);
            reader.setSkip(true);


            while (reader.hasNext() && count <= 600000) {

                try{
                    IAtomContainer molecule = reader.next();

                    molecule.setProperty("MOL_NUMBER_IN_FILE", Integer.toString(count));
                    molecule.setProperty("FILE_ORIGIN", file.getName().replace(".sdf", ""));

                    molecule.setProperty("SOURCE", source);




                    // Molecule original information

                    boolean foundOriginalSmiles = false;
                    molecule.setProperty("ORIGINAL_INCHI", "");
                    molecule.setProperty("ORIGINAL_INCHIKEY", "");

                    //trick to avoid having a molecule without even implicit hydrogens
                    try {
                            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                            CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
                            adder.addImplicitHydrogens(molecule);
                            System.out.println("Problem with adding implicit hydrogens and molecule configuration");
                            System.out.println(molecule);

                            Kekulization.kekulize(molecule);

                            System.out.println("Kekulization problem");
                            System.out.println(molecule);

                    }catch(CDKException e){
                        System.out.println("Problem with molecule in "+ source);
                        //System.out.println(molecule);
                    }

                    IAtomContainer tmpMolecule = molecule.clone();
                    for(Object p : molecule.getProperties().keySet()){

                        if(p.toString().toLowerCase().contains("smiles")){
                            tmpMolecule.setProperty("ORIGINAL_SMILES", molecule.getProperty(p));
                            foundOriginalSmiles = true;
                        }
                        if(p.toString().toLowerCase().contains("inchi") && !(p.toString().toLowerCase().contains("inchikey") || p.toString().toLowerCase().contains("inchi_key") || p.toString().toLowerCase().contains("inchi key")) ){
                            tmpMolecule.setProperty("ORIGINAL_INCHI", molecule.getProperty(p));
                        }
                        if(p.toString().toLowerCase().contains("inchikey") || p.toString().toLowerCase().contains("inchi_key") || p.toString().toLowerCase().contains("inchi key")){
                            tmpMolecule.setProperty("ORIGINAL_INCHIKEY", molecule.getProperty(p));
                        }
                    }

                    molecule = tmpMolecule;


                    if(!foundOriginalSmiles) {
                        molecule.setProperty("ORIGINAL_SMILES", smilesGenerator.create(molecule));
                    }




                    //Molecule curation
                    molecule = moleculeChecker.checkMolecule(molecule);


                    if (molecule != null) {

                        try {

                            List options = new ArrayList();
                            options.add(INCHI_OPTION.SNon);
                            options.add(INCHI_OPTION.ChiralFlagOFF);
                            options.add(INCHI_OPTION.AuxNone);
                            InChIGenerator gen = InChIGeneratorFactory.getInstance().getInChIGenerator(molecule, options );

                            molecule.setProperty("SIMPLE_INCHI", gen.getInchi());
                            molecule.setProperty("SIMPLE_INCHIKEY", gen.getInchiKey());
                        } catch (CDKException e) {

                            Integer totalBonds = molecule.getBondCount();
                            Integer ib = 0;
                            while (ib < totalBonds) {

                                IBond b = molecule.getBond(ib);
                                if (b.getOrder() == IBond.Order.UNSET) {
                                    //System.out.println(b.getOrder());
                                    b.setOrder(IBond.Order.SINGLE);

                                }
                                ib++;
                            }
                            List options = new ArrayList();
                            options.add(INCHI_OPTION.SNon);
                            options.add(INCHI_OPTION.ChiralFlagOFF);
                            options.add(INCHI_OPTION.AuxNone);
                            InChIGenerator gen = InChIGeneratorFactory.getInstance().getInChIGenerator(molecule, options );

                            molecule.setProperty("SIMPLE_INCHI", gen.getInchi());
                            molecule.setProperty("SIMPLE_INCHIKEY", gen.getInchiKey());
                        }


                        molecule.setProperty("SIMPLE_SMILES", smilesGenerator.create(molecule));


                        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd");
                        LocalDate localDate = LocalDate.now();

                        molecule.setProperty("ACQUISITION_DATE", dtf.format(localDate));

                        SourceNaturalProduct sourceNaturalProduct = ac2snp.createSNPlInstance(molecule);

                        sourceNaturalProduct.setContinent(databaseTypeChecker.checkContinent(this.source));

                        String taxa = databaseTypeChecker.checkKingdom(this.source);
                        if(taxa.equals("mixed")){
                            //do things db by db
                            if(source.equals("nubbedb")){
                                //there is a p at the beginning of each id for plants
                                if(molecule.getID().startsWith("p.")){
                                    taxa = "plants";
                                }else{
                                    taxa="animals";
                                }
                            }
                            else if(source.equals("npatlas")){
                                if(molecule.getID().startsWith("b")){
                                    taxa = "bacteria";
                                }else{
                                    taxa="fungi";
                                }
                            }
                            else{
                                taxa="notax";
                            }
                        }
                        sourceNaturalProduct.setOrganismText(new ArrayList<String>());
                        sourceNaturalProduct.organismText.add(taxa);


                        Hashtable<String, ArrayList<String>> sdfMetaData = searchMetaData(molecule);
                        if(sdfMetaData.containsKey("name")){
                            sourceNaturalProduct.setName(sdfMetaData.get("name").get(0));
                        }

                        if(sdfMetaData.containsKey("synonyms")){
                            sourceNaturalProduct.setSynonyms(sdfMetaData.get("synonyms"));
                        }

                        if(sdfMetaData.containsKey("citations")){
                            sourceNaturalProduct.setCitation(sdfMetaData.get("citations"));
                        }


                        if(!moleculeChecker.isForbiddenMolecule(molecule)){
                            sourceNaturalProductRepository.save(sourceNaturalProduct);
                        }
                    }



                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                count++;
                //System.out.println(count);

                if(count%50000==0){
                    System.out.println("Molecules read: "+count);
                }



            }




        } catch (IOException ex) {
            System.out.println("Oops ! File not found. Please check if the -in file or -out directory is correct");
            ex.printStackTrace();
            System.exit(0);
        }



    }

    @Override
    public ArrayList<IAtomContainer> returnCorrectMolecules() {
        return this.listOfMolecules;
    }

    @Override
    public String returnSource() {
        return this.source;
    }


    public Hashtable<String, ArrayList<String>> searchMetaData(IAtomContainer molecule){
        Hashtable<String, ArrayList<String>> foundMetaData = new Hashtable<>();


        for(Object p : molecule.getProperties().keySet()) {


            //search for name
            if (p.toString().toLowerCase().contains("name") && !(p.toString().toLowerCase().contains("database") )) {
                //check if list already created or not and add to it
                String n = molecule.getProperty(p);
                if(foundMetaData.containsKey("name")){
                    foundMetaData.get("name").add(n);
                }else{
                    foundMetaData.put("name" , new ArrayList<>());
                    foundMetaData.get("name").add(n);
                }

            }

            if(p.toString().toLowerCase().contains("synonym")){
                //see if need to split - split on ;

                //check if list already created or not and add to it
                String [] s = molecule.getProperty(p).toString().split(";");

                if(foundMetaData.containsKey("synonyms")){
                    for(String minis : s) {
                        foundMetaData.get("synonyms").add(minis);
                    }
                }else{
                    foundMetaData.put("synonyms" , new ArrayList<>());
                    for(String minis : s) {
                        foundMetaData.get("synonyms").add(minis);
                    }
                }

            }

            if(p.toString().toLowerCase().contains("pubmed_citation") || p.toString().toLowerCase().contains("citation") || p.toString().toLowerCase().contains("pubmed")
                || p.toString().toLowerCase().contains("doi") || p.toString().toLowerCase().contains("pmc")){

                // split on ;
                //check if list already created or not and add to it
                String [] cit = molecule.getProperty(p).toString().split(";");

                if(foundMetaData.containsKey("citations")){
                    for(String c : cit) {
                        foundMetaData.get("citations").add(c);
                    }
                }else{
                    foundMetaData.put("citations" , new ArrayList<>());
                    for(String c : cit) {
                        foundMetaData.get("citations").add(c);
                    }
                }
            }

        }




        return foundMetaData;
    }


}
