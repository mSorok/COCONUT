package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.function.IntBinaryOperator;

@Service
public class NPUnificationService {

    @Autowired
    SourceNaturalProductRepository sourceNaturalProductRepository;

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;


    PubchemFingerprinter pubchemFingerprinter = new PubchemFingerprinter( SilentChemObjectBuilder.getInstance() );

    CircularFingerprinter circularFingerprinter = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);


    SubstructureFingerprinter substructureFingerprinter = new SubstructureFingerprinter();

    ExtendedFingerprinter extendedFingerprinter = new ExtendedFingerprinter();



    private ArrayList<String> sourceNames;

    private Hashtable<String, String> sourceURLs;



    public void doWork(){

        System.out.println("NP unification InChi-key based");

        sourceNames = this.fetchSourceNames();

        sourceURLs = this.createSourceURLS();

        System.out.println("SOURCES  "+sourceNames);

        List<String> uniqueInchiKeys = sourceNaturalProductRepository.findUniqueInchiKeys();

        for(String oinchikey: uniqueInchiKeys){


            String inchikey =  oinchikey.split("\"")[3];  //oinchikey.toString();
            System.out.println(inchikey);
            List<SourceNaturalProduct> snpList = sourceNaturalProductRepository.findBySimpleInchiKey(inchikey);
            //System.out.println(snpList.get(0).simpleInchiKey+" "+snpList.get(0).source);

            //create a new UniqueNatural product

            UniqueNaturalProduct unp = new UniqueNaturalProduct();

            unp.setInchikey(inchikey);
            unp.setInchi(snpList.get(0).simpleInchi);
            unp.setSmiles(snpList.get(0).simpleSmiles);
            unp.setTotal_atom_number(snpList.get(0).getTotalAtomNumber());
            unp.setHeavy_atom_number(snpList.get(0).getHeavyAtomNumber());

            unp = uniqueNaturalProductRepository.save(unp);


            unp.name = "";

            unp.synonyms = new HashSet<>();
            unp.textTaxa = new HashSet<>();
            unp.taxid = new HashSet<>();
            unp.geoLocation = new HashSet<>();
            unp.citationDOI = new HashSet<>();
            unp.found_in_databases = new HashSet<>();
            unp.xrefs = new HashSet<>();
            unp.absolute_smiles = new Hashtable<>();

            //associate the UniqueNaturalProduct entry to each of the sources
            for(SourceNaturalProduct snp : snpList){
                snp.setUniqueNaturalProduct(unp);
                sourceNaturalProductRepository.save(snp);

                //add annotations from SourceNaturalProducts

                //name
                //checking if name doesn't contain DB name
               boolean nameIsReal = true;
                for(String dbname: this.sourceNames){
                    if(snp.getName() != null && (snp.getName().toLowerCase().contains(dbname) || snp.getName().startsWith("MLS") || snp.getName().startsWith("SMR") || snp.getName().contains("MLSMR"))){
                        nameIsReal=false;
                    }
                }


                if(nameIsReal && snp.getName() != null && ((unp.getName() == null || unp.getName() =="") &&  snp.getName().length()>3)){


                        String name = snp.getName().trim();

                        String[] names = name.split("\\\n");


                        unp.setName(names[0]);
                        if (names.length > 1) {
                            for (int i = 1; i < names.length; i++) {
                                unp.synonyms.add(names[i]);
                            }
                        }


                }
                else if( unp.getName() != null && unp.getName() != "" && snp.getName() != null){
                    if(snp.getSource().toLowerCase().contains("piellabdata")){
                        //replace name by ChebiName
                        unp.synonyms.add(unp.name);
                        unp.name = snp.getName().trim();

                        unp.nameTrustLevel=3;
                    }
                    else if(snp.getSource().toLowerCase().contains("chebi") && unp.nameTrustLevel<=2){

                        //replace name by ChebiName
                        unp.synonyms.add(unp.name);
                        unp.name = snp.getName().trim();

                        unp.nameTrustLevel=2;
                    }
                    else {
                        unp.synonyms.add(snp.getName().trim());
                    }
                }

                //synonyms
                if(snp.getSynonyms() != null){

                    String[] synonyms;

                    for(String sy : snp.getSynonyms()){

                        String[] names = sy.split("\\\n");
                        for (int i = 0; i < names.length; i++) {
                            unp.synonyms.add(names[i].trim());
                        }
                    }

                }


                //species
                if(snp.organismText != null ){

                    unp.textTaxa.addAll(snp.organismText);
                }
                if(snp.taxid != null){
                    unp.taxid.addAll(snp.taxid);
                }
                if(unp.textTaxa.size()>1 && unp.textTaxa.contains("notax")){
                    unp.textTaxa.remove("notax");
                }


                //geo
                if(snp.getGeographicLocation() != null){
                    unp.geoLocation.addAll(snp.getGeographicLocation());
                }
                if(snp.getContinent() != null){
                    unp.geoLocation.add(snp.getContinent());
                }
                if(unp.geoLocation.size()>1 && unp.geoLocation.contains("nogeo")){
                    unp.geoLocation.remove("nogeo");
                }


                //refs
                if(snp.getCitation() != null){
                    for(String cit : snp.getCitation()){
                        unp.citationDOI.add(cit.trim());
                    }

                }

                //cas
                if(snp.getCas() != null && snp.getCas() != ""){
                    unp.setCas(snp.getCas() );
                }

                //database
                if(snp.getSource() != null){
                    unp.found_in_databases.add(snp.getSource());

                    ArrayList<String> miniXref = new ArrayList<String>();
                    miniXref.add(snp.getSource());
                    miniXref.add(snp.idInSource);
                    if(sourceURLs.containsKey(snp.getSource())) {
                        miniXref.add(sourceURLs.get(snp.getSource()));
                    }
                    unp.xrefs.add(miniXref);
                }

                //Absolute smiles (with stereochemistry)
                if(snp.getAbsoluteSmiles() != null && !snp.getAbsoluteSmiles().equals("")) {
                    if (unp.absolute_smiles.containsKey(snp.getAbsoluteSmiles())) {
                        unp.absolute_smiles.get(snp.getAbsoluteSmiles()).add(snp.getSource());

                    } else {
                        HashSet newSourceList = new HashSet();
                        newSourceList.add(snp.getSource());
                        unp.absolute_smiles.put(snp.getAbsoluteSmiles(), newSourceList);
                    }
                }
                else{
                    if(unp.absolute_smiles.containsKey("nostereo")) {
                        unp.absolute_smiles.get("nostereo").add(snp.getSource());
                    }
                    else{
                        HashSet newSourceList = new HashSet();
                        newSourceList.add(snp.getSource());
                        unp.absolute_smiles.put("nostereo", newSourceList);
                    }
                }

            }

            unp = uniqueNaturalProductRepository.save(unp);

            //compute molecular parameters for the UniqueNaturalProduct
            unp = computeFingerprints(unp);
            unp = computeAdditionalMolecularFeatures(unp);
            uniqueNaturalProductRepository.save(unp);

        }



    }









    public UniqueNaturalProduct computeAdditionalMolecularFeatures(UniqueNaturalProduct m){
        AllRingsFinder arf = new AllRingsFinder();
        MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
        AtomContainerManipulator acm = new AtomContainerManipulator();

        IAtomContainer im = atomContainerToUniqueNaturalProductService.createAtomContainer(m);


        // count rings
        try {
            IRingSet rs = arf.findAllRings(im, 20);

            m.setMax_number_of_rings(rs.getAtomContainerCount());

            Cycles   cycles = Cycles.sssr(im);
            IRingSet rings  = cycles.toRingSet();
            m.setMin_number_of_rings(rings.getAtomContainerCount()); //SSSR


        } catch (CDKException e) {
            System.out.println("Too complex: "+m.getSmiles());
        }

        //compute molecular formula
        m.setMolecular_formula(mfm.getString(mfm.getMolecularFormula(im) ));

        //compute number of carbons, of nitrogens, of oxygens
        m.setNumber_of_carbons(mfm.getElementCount(mfm.getMolecularFormula(im), "C"));

        m.setNumber_of_oxygens(mfm.getElementCount(mfm.getMolecularFormula(im), "O"));

        m.setNumber_of_nitrogens(mfm.getElementCount(mfm.getMolecularFormula(im), "N"));

        m.setMolecular_weight( acm.getMass(im) );


        // cleaning the NaNs
        if(m.getMolecular_weight().isNaN()){
            m.setMolecular_weight(0.0);
        }


        //get bond count
        IBond[] bonds = acm.getBondArray(im);
        int bondCount = 0;
        for(IBond b : bonds){
            if(b.getAtomCount() == 2){
                if(!b.getAtom(0).getSymbol().equals("H") && !b.getAtom(1).getSymbol().equals("H")){
                    bondCount++;
                }
            }
        }
        m.setBond_count(bondCount);




        return(m);
    }




    public ArrayList fetchSourceNames(){
        List<Object> uniqueSourceNames = sourceNaturalProductRepository.findUniqueSourceNames();
        ArrayList<String> tmpArray = new ArrayList<>();

        for(Object usn: uniqueSourceNames) {
            tmpArray.add(usn.toString().toLowerCase());
        }
        return tmpArray;
    }





    public UniqueNaturalProduct computeFingerprints(UniqueNaturalProduct np){



        IAtomContainer ac = atomContainerToUniqueNaturalProductService.createAtomContainer(np);

        // Addition of implicit hydrogens & atom typer
        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(ac.getBuilder());
        CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(ac.getBuilder() );
        for (int j = 0; j < ac.getAtomCount(); j++) {
            IAtom atom = ac.getAtom(j);
            IAtomType type = null;
            try {
                type = matcher.findMatchingAtomType(ac, atom);
                AtomTypeManipulator.configure(atom, type);
            } catch (CDKException e) {
                e.printStackTrace();
            }

        }

        try {
            adder.addImplicitHydrogens(ac);
        } catch (CDKException e) {
            e.printStackTrace();
        }
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
        AtomContainerManipulator.removeNonChiralHydrogens(ac);

        try {

            String s = pubchemFingerprinter.getBitFingerprint(ac).asBitSet().toString();
            ArrayList<Integer> pcl = new ArrayList<>();
            s = s.replace(" ", "");s = s.replace("\"", "");s = s.replace("{", "");s = s.replace("}", "");
            String [] sl = s.split(",");
            for(String c : sl){
                try {
                    pcl.add(Integer.parseInt(c));
                }catch (NumberFormatException e){ e.printStackTrace(); }
            }
            np.setPubchemFingerprint(pcl);

            np.pubfp = new HashMap<>();
            np.pubfp.put(new Integer(pcl.size()), pcl);


            s = circularFingerprinter.getBitFingerprint(ac).asBitSet().toString();
            pcl = new ArrayList<>();
            s = s.replace(" ", "");s = s.replace("\"", "");s = s.replace("{", "");s = s.replace("}", "");
            sl = s.split(",");
            for(String c : sl){
                try {
                    pcl.add(Integer.parseInt(c));
                }catch (NumberFormatException e){ e.printStackTrace(); }
            }
            np.setCircularFingerprint(pcl);

            s = extendedFingerprinter.getBitFingerprint(ac).asBitSet().toString();
            pcl = new ArrayList<>();
            s = s.replace(" ", "");s = s.replace("\"", "");s = s.replace("{", "");s = s.replace("}", "");
            sl = s.split(",");
            for(String c : sl){
                try {
                    pcl.add(Integer.parseInt(c));
                }catch (NumberFormatException e){ e.printStackTrace(); }
            }
            np.setExtendedFingerprint(pcl);



            //Bits and String for PubChem

            try {
                //for PubChem

                BitSet bitsOn = pubchemFingerprinter.getBitFingerprint(ac).asBitSet();
                String pubchemBitString = "";

                for (int i = 0; i <= bitsOn.length(); i++) {
                    if (bitsOn.get(i)) {
                        pubchemBitString += "1";
                    } else {
                        pubchemBitString += "0";
                    }
                }

                np.setPubchemBits(bitsOn.toByteArray());
                np.setPubchemBitsString(pubchemBitString);



            } catch (CDKException | UnsupportedOperationException e) {
                e.printStackTrace();
            }


        } catch (CDKException | UnsupportedOperationException e) {
            e.printStackTrace();
        }
        return np;
    }


    public Hashtable createSourceURLS(){
        Hashtable<String, String> urls = new Hashtable<>();

        if( this.sourceNames != null && !this.sourceNames.isEmpty()){
            for(String sourceName : this.sourceNames){
                /*if(sourceName.equals("biofacquim")){
                    urls.put("biofacquim", "https://biofacquim.herokuapp.com/");
                }*/

                if(sourceName.equals("bitterdb")){
                    urls.put("bitterdb", "http://bitterdb.agri.huji.ac.il/compound.php?id=");
                }

                if(sourceName.equals("carotenoids")){
                    urls.put("carotenoids", "http://carotenoiddb.jp/Entries/");
                }

                if(sourceName.equals("chebi_np")){
                    urls.put("chebi_np", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:");
                }

                if(sourceName.equals("chembl_np")){
                    urls.put("chembl_np", "https://www.ebi.ac.uk/chembl/compound_report_card/");
                }


                if(sourceName.equals("cmaup")){
                    urls.put("cmaup", "http://bidd2.nus.edu.sg/CMAUP/searchresults.php?keyword_search=");
                }

                if(sourceName.equals("pubchem_tested_np")){
                    urls.put("pubchem_tested_np", "https://pubchem.ncbi.nlm.nih.gov/compound/");
                }

                if(sourceName.equals("drugbanknp")){
                    urls.put("drugbanknp", "https://www.drugbank.ca/drugs/");
                }


                if(sourceName.equals("chemspidernp")){
                    urls.put("chemspidernp", "http://www.chemspider.com/Chemical-Structure.");
                }

                if(sourceName.equals("np_atlas_2019_12") ){
                    urls.put("np_atlas_2019_12", "https://www.npatlas.org/joomla/index.php/explore/compounds#npaid=");
                }

                if(sourceName.equals("npatlas")){
                    urls.put("npatlas", "https://www.npatlas.org/joomla/index.php/explore/compounds#npaid=");
                }




                if(sourceName.equals("exposome-explorer")){
                    urls.put("exposome-explorer", "http://exposome-explorer.iarc.fr/compounds/");
                }

                if(sourceName.equals("fooddb")){
                    urls.put("fooddb", "https://foodb.ca/compounds/");
                }

                if(sourceName.equals("knapsack")){
                    urls.put("knapsack", "http://www.knapsackfamily.com/knapsack_core/information.php?mode=r&word=");
                }

                if(sourceName.equals("npass")){
                    urls.put("npass", "http://bidd2.nus.edu.sg/NPASS/browse_np.php?compound=");
                }

                if(sourceName.equals("nubbe")){
                    urls.put("nubbe", "https://nubbe.iq.unesp.br/portal/nubbe-search.html");
                }

                if(sourceName.equals("phenolexplorer")){
                    urls.put("phenolexplorer", "http://phenol-explorer.eu/compounds/");
                }


                if(sourceName.equals("sancdb")){
                    urls.put("sancdb", "https://sancdb.rubi.ru.ac.za/compounds/");
                }


                if(sourceName.equals("supernatural2")){
                    urls.put("supernatural2", "http://bioinf-applied.charite.de/supernatural_new/index.php?site=compound_search&id=");
                }


                if(sourceName.equals("tcmdb_taiwan")){
                    urls.put("tcmdb_taiwan", "http://tcm.cmu.edu.tw/");
                }

                if(sourceName.equals("tppt")){
                    urls.put("tppt", "https://www.agroscope.admin.ch/agroscope/en/home/publications/apps/tppt.html");
                }

                if(sourceName.equals("vietherb")){
                    urls.put("vietherb", "https://vietherb.com.vn/metabolites/");
                }

                if(sourceName.equals("streptomedb")){
                    urls.put("streptomedb", "http://132.230.56.4/streptomedb2/get_drugcard/");
                }






            }
        }



        return urls;

    }



}
