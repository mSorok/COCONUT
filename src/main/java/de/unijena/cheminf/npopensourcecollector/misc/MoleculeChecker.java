package de.unijena.cheminf.npopensourcecollector.misc;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.surface.NeighborList;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.springframework.stereotype.Service;

import javax.sound.midi.SysexMessage;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

@Service
public class MoleculeChecker {

    private static final int MIN_HEAVY_ATOM_COUNT = 5;
    private static final int MAX_HEAVY_ATOM_COUNT = 210;



    private final String[] check = {"C", "H", "N", "O", "P", "S", "Cl", "F", "As", "Se", "Br", "I", "B", "Na", "Si", "K", "Fe"};
    private final HashSet<String> symbols2Check = new HashSet<String>(Arrays.asList(check));

    private final String[] forbiddenInchiKeys = {"OOHPORRAEMMMCX-UHFFFAOYSA-N", "ATSPGPYEGAPMOB-UHFFFAOYSA-N"};
    private final  HashSet<String> inchis2Check = new HashSet<String>(Arrays.asList(forbiddenInchiKeys));



    MoleculeConnectivityChecker mcc;



    public IAtomContainer checkMolecule(IAtomContainer molecule){

        IAtomContainer oriMol = molecule;


        mcc = BeanUtil.getBean(MoleculeConnectivityChecker.class);

        SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Absolute);
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        if(!containsStrangeElements(molecule)) {

            /**
             * Checking for connectivity and selecting the biggest component
             */

            List<IAtomContainer> listAC = mcc.checkConnectivity(molecule);
            if( listAC.size()>=1 ){
                IAtomContainer biggestComponent = listAC.get(0);
                for(IAtomContainer partac : listAC){
                    if(partac.getAtomCount()>biggestComponent.getAtomCount()){
                        biggestComponent = partac;
                    }
                }
                molecule = biggestComponent;

                int nbheavyatoms = 0;
                for(IAtom a : molecule.atoms()){
                    if(!a.getSymbol().equals("H")){
                        nbheavyatoms++;
                    }
                }
                if(nbheavyatoms<= MIN_HEAVY_ATOM_COUNT || nbheavyatoms>=MAX_HEAVY_ATOM_COUNT){
                    return null;
                }
            }
            else{

                return null;
            }



            // check ID

            if (molecule.getID() == "" || molecule.getID() == null) {
                for (Object p : molecule.getProperties().keySet()) {

                    if (p.toString().toLowerCase().contains("id")) {
                        molecule.setID(molecule.getProperty(p.toString()));

                    }

                }
                if (molecule.getID() == "" || molecule.getID() == null) {
                    molecule.setID(molecule.getProperty("MOL_NUMBER_IN_FILE"));
                    //this.molecule.setProperty("ID", this.molecule.getProperty("MOL_NUMBER_IN_FILE"));
                }


            }

            Map<Object, Object> properties = molecule.getProperties();
            String id = molecule.getID();
            oriMol = molecule;


            //System.out.println(id);
            molecule = sanitizeMolecularBonds(molecule);

            if(molecule == null){
                //System.out.println("No sanitization possible");
                return null;
            }
/*
            try {
                System.out.println(sg.create(molecule));
            } catch (CDKException e) {
                e.printStackTrace();
            }
*/
            //Normalizing the ionization states

            //python3 GetParentSourceNP.py "C[C@H](O)c1ccccc1"

            try {
                //String command = "evaluate -e majorMicrospecies('7.4') "+sg.create(molecule);
                String command = "python3 GetParentSourceNP.py "+sg.create(molecule);
                Process process = Runtime.getRuntime().exec(command);

                BufferedReader reader = new BufferedReader(
                        new InputStreamReader(process.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    String correctSmi = line.replace("\n", "");
                    if (!line.contains("Normalizer") && !line.contains("charge") && !line.equals("") ){

                        System.out.println(correctSmi);

                    try {
                        molecule = sp.parseSmiles(correctSmi);
                        molecule.setProperties(properties);
                        molecule.setID(id);
                    } catch (CDKException | IllegalArgumentException e) {
                        e.printStackTrace();

                        molecule = oriMol;
                    }
                }
                }

                reader.close();


            }catch(Exception e){
                System.out.println("Couldn't start Major Microscpecies calculation: something went wrong");
                e.printStackTrace();
                molecule = oriMol;
            }






            //ElectronDonation model = ElectronDonation.cdk();
            //CycleFinder cycles = Cycles.cdkAromaticSet();
            //Aromaticity aromaticity = new Aromaticity(model, cycles);





            try {
                String smi = sg.create(molecule);
                molecule = sp.parseSmiles(smi);
                molecule.setProperties(properties);
                molecule.setID(id);
            } catch (CDKException | IllegalArgumentException e) {
                e.printStackTrace();
            }


            // Addition of implicit hydrogens & atom typer
            CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(molecule.getBuilder());
            for (int j = 0; j < molecule.getAtomCount(); j++) {
                IAtom atom = molecule.getAtom(j);
                IAtomType type = null;
                try {
                    type = matcher.findMatchingAtomType(molecule, atom);
                    AtomTypeManipulator.configure(atom, type);
                } catch (CDKException e) {
                    e.printStackTrace();
                }

            }


            CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder() );

            try {
                adder.addImplicitHydrogens(molecule);
            } catch (CDKException e) {
                e.printStackTrace();
            }

            AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
            AtomContainerManipulator.removeNonChiralHydrogens(molecule);




            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(molecule);
                //Adding aromaticity to molecules when needed
                //aromaticity.apply(molecule);

            } catch (CDKException e) {
                e.printStackTrace();
            }




            //Fixing molecular bonds
            try {
                Kekulization.kekulize(molecule);

            } catch (CDKException e1) {
                //e1.printStackTrace();
            } catch (IllegalArgumentException e) {
                //System.out.println("Could not kekulize molecule "+ this.molecule.getID());
            }





            return molecule;
        }
        return null;
    }



    private boolean containsStrangeElements(IAtomContainer molecule) {
        if(molecule.getAtomCount()>0) {
            for (IAtom atom : molecule.atoms()) {
                if (!symbols2Check.contains(atom.getSymbol())) {
                    System.out.println("contains strange");
                    System.out.println(atom.getSymbol());
                    System.out.println(molecule.getID());
                    return true;
                }
            }
        }
        return false;
    }

    public boolean isForbiddenMolecule(IAtomContainer molecule){
        String inchikey = molecule.getProperty("INCHIKEY");
        if(inchis2Check.contains(inchikey)){
            return true;
        }
        return false;

    }

    public IAtomContainer removeExcessiveDoubleBonds(IAtomContainer atomContainer, IAtom atom, Integer maxValency){

        List<IBond> listOfBonds = atomContainer.getConnectedBondsList(atom);


        int b = 0;
        while(AtomContainerManipulator.getBondOrderSum(atomContainer, atom)>maxValency && b<listOfBonds.size()){
            listOfBonds.get(b).setOrder(IBond.Order.SINGLE);
            b++;
        }

        return atomContainer;

    }


    public IAtomContainer removeExtraHydrogen(IAtomContainer acToModify, IAtom atomToModify, Integer maxValency){




        List<IAtom> neighborsOfTheSuspiciousAtom = acToModify.getConnectedAtomsList(atomToModify);



        if(atomToModify.getSymbol().equals("C")){
            // check if no N or O or P in the neighbors. If yes, return, we'll be removing bonds and hydrogens from them instead
            int t = 0;
            while(t<neighborsOfTheSuspiciousAtom.size()){
                if(neighborsOfTheSuspiciousAtom.get(t).getSymbol().equals("O") || neighborsOfTheSuspiciousAtom.get(t).getSymbol().equals("N") || neighborsOfTheSuspiciousAtom.get(t).getSymbol().equals("P") ){
                    return acToModify;
                }
                t++;
            }

        }

        acToModify = removeExcessiveDoubleBonds(acToModify, atomToModify, maxValency);


        int nbHydrogens = atomToModify.getImplicitHydrogenCount();


        while(AtomContainerManipulator.getBondOrderSum(acToModify, atomToModify)>maxValency && nbHydrogens>0){

            nbHydrogens = nbHydrogens-1;
            atomToModify.setImplicitHydrogenCount(nbHydrogens);

        }

        if(AtomContainerManipulator.getBondOrderSum(acToModify, atomToModify)>maxValency){
            acToModify = null;
            return acToModify;
        }

        return(acToModify);
    }


    IAtomContainer sanitizeMolecularBonds(IAtomContainer atomContainer){



        // Run sanitization on the number of bonds expected for carbons, nitrogens and oxygens
        // and
        //Homogenize pseudo atoms - all pseudo atoms (PA) as a "*"
        for (int u = 1; u < atomContainer.getAtomCount(); u++) {



            if (atomContainer.getAtom(u) instanceof IPseudoAtom) {

                atomContainer.getAtom(u).setSymbol("*");
                atomContainer.getAtom(u).setAtomTypeName("X");
                ((IPseudoAtom) atomContainer.getAtom(u)).setLabel("*");

            }

            //check if the number of bonds is correct. If it's too big, remove the attached hydrogens
            double bondCount = AtomContainerManipulator.getBondOrderSum(atomContainer, atomContainer.getAtom(u));


            if(atomContainer.getAtom(u).getSymbol().equals("Br") && bondCount>1 ){

                atomContainer = removeExtraHydrogen(atomContainer, atomContainer.getAtom(u), 1);
            }
            else if(atomContainer.getAtom(u).getSymbol().equals("O") && bondCount>2 ){
                atomContainer = removeExtraHydrogen(atomContainer, atomContainer.getAtom(u), 2);
            }
            else if(atomContainer.getAtom(u).getSymbol().equals("N") && bondCount>3 ){
                atomContainer = removeExtraHydrogen(atomContainer, atomContainer.getAtom(u), 3);

            }
            else if(atomContainer.getAtom(u).getSymbol().equals("C") && bondCount>4 ){
                atomContainer = removeExtraHydrogen(atomContainer, atomContainer.getAtom(u), 4);
            }
            else if(atomContainer.getAtom(u).getSymbol().equals("P") && bondCount>5 ){
                atomContainer = removeExtraHydrogen(atomContainer, atomContainer.getAtom(u), 5);
            }

            if(atomContainer == null){
                return null;
            }

        }





        return atomContainer;
    }


}
