package de.unijena.cheminf.npopensourcecollector.readers;

import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.DatabaseTypeChecker;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeChecker;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.services.AtomContainerToSourceNaturalProductService;
import net.sf.jniinchi.INCHI_OPTION;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;

public class MOLReader implements Reader {

    File file;
    ArrayList<IAtomContainer> listOfMolecules;

    private IteratingSDFReader reader = null;

    SourceNaturalProductRepository sourceNaturalProductRepository;
    AtomContainerToSourceNaturalProductService ac2snp;


    MoleculeChecker moleculeChecker;
    DatabaseTypeChecker databaseTypeChecker;

    String source;

    public MOLReader(){
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

        this.source = file.getName().toLowerCase().replace(".mol", "");



        try{

            reader = new IteratingSDFReader(new FileInputStream(file), DefaultChemObjectBuilder.getInstance());
            System.out.println("MOL reader creation and inserting in MongoDB for "+source);
            reader.setSkip(true);


            while (reader.hasNext() && count <= 600000) {

                try{
                    IAtomContainer molecule = reader.next();

                    molecule.setProperty("MOL_NUMBER_IN_FILE", Integer.toString(count));
                    molecule.setProperty("FILE_ORIGIN", file.getName().replace(".mol", ""));

                    molecule.setProperty("SOURCE", source);



                    // Molecule original information

                    boolean foundOriginalSmiles = false;
                    molecule.setProperty("ORIGINAL_INCHI", "");
                    molecule.setProperty("ORIGINAL_INCHIKEY", "");

                    //trick to avoid having a molecule without even implicit hydrogens
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                    CDKHydrogenAdder adder =
                            CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
                    adder.addImplicitHydrogens(molecule);

                    IAtomContainer tmpMolecule = molecule.clone();
                    for(Object p : molecule.getProperties().keySet()){

                        if(p.toString().toLowerCase().contains("smiles")){
                            tmpMolecule.setProperty("ORIGINAL_SMILES", molecule.getProperty(p));
                            foundOriginalSmiles = true;
                        }
                        if(p.toString().toLowerCase().contains("inchi") && !(p.toString().toLowerCase().contains("inchikey") || p.toString().toLowerCase().contains("inchi_key")) ){
                            tmpMolecule.setProperty("ORIGINAL_INCHI", molecule.getProperty(p));
                        }
                        if(p.toString().toLowerCase().contains("inchikey") || p.toString().toLowerCase().contains("inchi_key")){
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


                        if(!moleculeChecker.isForbiddenMolecule(molecule)){
                            sourceNaturalProductRepository.save(sourceNaturalProduct);
                        }
                    }



                } catch (Exception ex) {
                    //ex.printStackTrace();
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
}
