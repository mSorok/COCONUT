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
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SMILESReader implements Reader {


    File file;
    ArrayList<IAtomContainer> listOfMolecules;
    private LineNumberReader smilesReader;
    SourceNaturalProductRepository sourceNaturalProductRepository;
    AtomContainerToSourceNaturalProductService ac2snp;
    MoleculeChecker moleculeChecker;
    DatabaseTypeChecker databaseTypeChecker;
    String source;


    public SMILESReader(){
        this.listOfMolecules = new ArrayList<IAtomContainer>();
        sourceNaturalProductRepository = BeanUtil.getBean(SourceNaturalProductRepository.class);
        ac2snp = BeanUtil.getBean(AtomContainerToSourceNaturalProductService.class);
        moleculeChecker = BeanUtil.getBean(MoleculeChecker.class);
        databaseTypeChecker = BeanUtil.getBean(DatabaseTypeChecker.class);

    }




    @Override
    public void readFile(File file) {
        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique );


        this.file = file;
        int count = 1;
        String line;

        this.source = file.getName().toLowerCase().replace(".smi", "");

        try(BufferedReader smilesReader = new BufferedReader(new FileReader(file))) {
            System.out.println("SMILES reader creation and inserting in MongoDB for "+source);


            while ((line = smilesReader.readLine()) != null  && count <= 600000) {

                String smiles_names = line;

                if(!line.contains("smiles")) {
                    try {
                        String[] splitted = smiles_names.split("\\s+"); //splitting the canonical smiles format: SMILES \s mol name
                        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

                        IAtomContainer molecule = null;
                        try {
                            molecule = sp.parseSmiles(splitted[0]);


                            molecule.setProperty("MOL_NUMBER_IN_FILE", Integer.toString(count));
                            if(splitted.length==2) {
                                molecule.setProperty("ID", splitted[1]);
                                molecule.setID(splitted[1]);
                            }else if(splitted.length==1){
                                //no id
                                molecule.setProperty("ID", Integer.toString(count));
                                molecule.setID(Integer.toString(count));
                            }else if(splitted.length>2){
                                //join everything after 1
                                ArrayList<String> splitted2 = new ArrayList<String>(Arrays.asList(splitted));

                                String nid = String.join(" ", splitted2);
                                molecule.setProperty("ID", nid);
                                molecule.setID(nid);
                            }

                            molecule.setProperty("FILE_ORIGIN", file.getName().replace(".smi", ""));

                            molecule.setProperty("SOURCE", source);

                            molecule.setProperty("ORIGINAL_INCHI", "");
                            molecule.setProperty("ORIGINAL_INCHIKEY", "");
                            molecule.setProperty("ORIGINAL_SMILES", splitted[0]);


                            molecule = moleculeChecker.checkMolecule(molecule);

                            if (molecule != null){
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

                        } catch (InvalidSmilesException e) {
                            e.printStackTrace();
                            System.out.println(line);
                            System.out.println(splitted);
                            System.out.println(splitted[0]);
                            smilesReader.skip(count - 1);
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    count++;

                    if(count%50000==0){
                        System.out.println("Molecules read: "+count);
                    }
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
