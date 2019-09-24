package de.unijena.cheminf.npopensourcecollector.readers;

import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.DatabaseTypeChecker;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeChecker;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.services.AtomContainerToSourceNaturalProductService;
import net.sf.jniinchi.INCHI_OPTION;
import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
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

public class CSVReader implements Reader {


    File file;
    ArrayList<IAtomContainer> listOfMolecules;
    private LineNumberReader inchiReader;
    SourceNaturalProductRepository sourceNaturalProductRepository;
    AtomContainerToSourceNaturalProductService ac2snp;
    MoleculeChecker moleculeChecker;
    DatabaseTypeChecker databaseTypeChecker;
    String source;

    public CSVReader(){
        this.listOfMolecules = new ArrayList<IAtomContainer>();
        sourceNaturalProductRepository = BeanUtil.getBean(SourceNaturalProductRepository.class);
        ac2snp = BeanUtil.getBean(AtomContainerToSourceNaturalProductService.class);
        moleculeChecker = BeanUtil.getBean(MoleculeChecker.class);
        databaseTypeChecker = BeanUtil.getBean(DatabaseTypeChecker.class);
    }


    @Override
    public void readFile(File file) {

        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique );
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        this.file = file;
        this.source = file.getName().toLowerCase().replace(".csv", "");


        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(this.file));
            //read the header

            // if the first line is the header
            ArrayList<String> header = new ArrayList<String>(Arrays.asList( bufferedReader.readLine().split(",")));
            //System.out.println(header);

            Integer indexOfID=null;
            Integer indexOfName=null;
            Integer indexOfReference=null;
            Integer indexOfCitation=null;
            Integer indexOfDOI=null;
            Integer indexOfSMILES=null;
            Integer indexOfInchi=null;
            Integer indexOfInchikey=null;
            Integer indexOfKingdom=null;
            Integer indexOfGenus=null;
            Integer indexOfSpecies=null;
            Integer indexOfGeo=null;
            Integer indexOfCode=null;

            for(String item : header){

                if( item.toLowerCase().equals("id") || item.toLowerCase().equals("identifier") ){
                    indexOfID = header.indexOf(item);
                }
                if(item.toLowerCase().contains("name")){
                    indexOfName = header.indexOf(item);
                }
                if(item.toLowerCase().contains("ref")){
                    indexOfReference = header.indexOf(item);
                }
                if(item.toLowerCase().contains("citation")){
                    indexOfCitation = header.indexOf(item);
                }
                if(item.toLowerCase().contains("doi")){
                    indexOfDOI = header.indexOf(item);
                }
                if(item.toLowerCase().contains("smiles")){
                    indexOfSMILES = header.indexOf(item);
                }
                if(item.toLowerCase().contains("inchi") && !item.toLowerCase().contains("inchikey")){
                    indexOfInchi = header.indexOf(item);
                }
                if(item.toLowerCase().contains("inchikey")){
                    indexOfInchikey = header.indexOf(item);
                }
                if(item.toLowerCase().contains("kingdom")){
                    indexOfKingdom = header.indexOf(item);
                }
                if(item.toLowerCase().contains("genu")){
                    indexOfGenus = header.indexOf(item);
                }
                if(item.toLowerCase().contains("specie")){
                    indexOfSpecies = header.indexOf(item);
                }
                if(item.toLowerCase().contains("geo") || item.toLowerCase().contains("site") || item.toLowerCase().contains("local")){
                    indexOfGeo = header.indexOf(item);
                }
                if( item.toLowerCase().contains("code") || item.toLowerCase().contains(this.source) ){
                    indexOfCode = header.indexOf(item);
                }

            }


            if(indexOfID== null && indexOfCode != null){
                indexOfID = indexOfCode;
            }


            //read the rest of the file
            int count = 1;
            String line;

            while ((line = bufferedReader.readLine()) != null  && count <= 600000) {

                ArrayList<String> dataline = new ArrayList<String>(Arrays.asList( line.split(",")));
                try{

                    IAtomContainer molecule = null;

                    if(indexOfSMILES != null){
                        molecule = sp.parseSmiles(dataline.get(indexOfSMILES));

                        molecule.setProperty("FILE_ORIGIN", file.getName().replace(".csv", ""));
                        molecule.setProperty("SOURCE", source);
                        molecule.setProperty("ORIGINAL_SMILES", dataline.get(indexOfSMILES));


                        if(indexOfInchi != null){
                            molecule.setProperty("ORIGINAL_INCHI", dataline.get(indexOfInchi));

                        }
                        if(indexOfInchikey != null){
                            molecule.setProperty("ORIGINAL_INCHIKEY", dataline.get(indexOfInchikey));
                        }

                    }else if(indexOfInchi != null){
                        // READING InCHI
                        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                        InChIToStructure intostruct = factory.getInChIToStructure(dataline.get(indexOfInchi), DefaultChemObjectBuilder.getInstance());
                        INCHI_RET ret = intostruct.getReturnStatus();
                        if (ret == INCHI_RET.WARNING) {
                            // Structure generated, but with warning message
                            System.out.println("InChI warning: " + intostruct.getMessage());
                        } else if (ret != INCHI_RET.OKAY) {
                            // Structure generation failed
                            throw new CDKException("Structure generation failed failed: " + ret.toString() + " [" + intostruct.getMessage() + "]");
                        }

                        molecule = intostruct.getAtomContainer();
                        if(indexOfInchikey != null){
                            molecule.setProperty("ORIGINAL_INCHIKEY", dataline.get(indexOfInchikey));
                        }
                    }

                    if(molecule != null){
                        if(indexOfID != null){
                            molecule.setID(dataline.get(indexOfID));
                            molecule.setProperty("ID", dataline.get(indexOfID));
                        }
                        else if(indexOfName != null){
                            molecule.setID(dataline.get(indexOfName));
                            molecule.setProperty("ID", dataline.get(indexOfName));
                        }
                        else{
                            molecule.setID(Integer.toString(count));
                            molecule.setProperty("ID", Integer.toString(count));
                        }

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
                                else if(source.equals("biofacquim")){
                                    taxa = dataline.get(indexOfKingdom);
                                }
                                else{
                                    taxa="notax";
                                }
                            }
                            sourceNaturalProduct.setOrganismText(new ArrayList<String>());
                            sourceNaturalProduct.organismText.add(taxa);

                            if(indexOfKingdom != null){
                                sourceNaturalProduct.organismText.add(dataline.get(indexOfKingdom));
                            }
                            if(indexOfGenus != null){
                                sourceNaturalProduct.organismText.add(dataline.get(indexOfGenus));
                            }
                            if(indexOfSpecies != null){
                                sourceNaturalProduct.organismText.add(dataline.get(indexOfSpecies));
                            }

                            //GEOGRAPHY
                            sourceNaturalProduct.setContinent(databaseTypeChecker.checkContinent(this.source));
                            if(indexOfGeo != null){
                                sourceNaturalProduct.geographicLocation = new ArrayList<>();
                                sourceNaturalProduct.geographicLocation.add(dataline.get(indexOfGeo));
                            }

                            //citation reference and doi
                            if(indexOfCitation != null || indexOfDOI != null || indexOfReference != null){
                                sourceNaturalProduct.citation = new ArrayList<>();
                                if(indexOfCitation != null){
                                    sourceNaturalProduct.citation.add(dataline.get(indexOfCitation));
                                }
                                if(indexOfDOI != null){
                                    sourceNaturalProduct.citation.add(dataline.get(indexOfDOI));
                                }
                                if(indexOfReference != null){
                                    sourceNaturalProduct.citation.add(dataline.get(indexOfReference));
                                }
                            }

                            if(!moleculeChecker.isForbiddenMolecule(molecule)){
                                sourceNaturalProductRepository.save(sourceNaturalProduct);
                            }
                        }
                    }
                    else{
                        System.out.println("No molecular structure detected");
                    }




                } catch (CDKException e ) {
                    e.printStackTrace();
                    System.out.println(line);
                }

                count++;
            }


        } catch (IOException e ) {
            e.printStackTrace();
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
