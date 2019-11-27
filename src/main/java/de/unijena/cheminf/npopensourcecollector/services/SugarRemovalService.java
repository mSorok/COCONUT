package de.unijena.cheminf.npopensourcecollector.services;


import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeConnectivityChecker;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Ring;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.springframework.stereotype.Service;

import java.util.*;



@Service
public class SugarRemovalService {



    public List<IAtomContainer> linearSugars;
    public List<IAtomContainer> ringSugars;
    public List<DfPattern> patternListLinearSugars;


    UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();
    MoleculeConnectivityChecker mcc;





    public void getSugarChains() {

        linearSugars = new ArrayList<IAtomContainer>();
        ringSugars = new ArrayList<IAtomContainer>();
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String[] linearSugarsSmilesList = {"C(C(C(C(C(C=O)O)O)O)O)O", "C(C(CC(C(CO)O)O)O)(O)=O", "C(C(C(CC(=O)O)O)O)O",
                "C(C(C(C(C(CO)O)O)O)=O)O", "C(C(C(C(C(CO)O)O)O)O)O", "C(C(C(C(CC=O)O)O)O)O",
                "OCC(O)C(O)C(O)C(O)CO", "O=CC(O)C(O)C(O)C(O)CO",  "CCCCC(O)C(=O)O", "CC(=O)CC(=O)CCC(=O)O",
                "O=C(O)CC(O)CC(=O)O", "O=C(O)C(=O)C(=O)C(O)C(O)CO",
                "O=C(O)CCC(O)C(=O)O", "O=CC(O)C(O)C(O)C(O)CO", "O=C(CO)C(O)C(O)CO"};



        String [] ringSugarsSmilesList = { "C1CCOC1", "C1CCOCC1", "C1CCCOCC1"};



        //adding linear sugars
        for (String smiles : linearSugarsSmilesList) {
            try {
                linearSugars.add(smilesParser.parseSmiles(smiles));
            } catch (InvalidSmilesException ex) {
                System.out.println("could not read linear sugar");
            }
        }


        //adding ring sugars
        for (String smiles : ringSugarsSmilesList) {
            try {
                ringSugars.add(smilesParser.parseSmiles(smiles));
            } catch (InvalidSmilesException ex) {
                System.out.println("could not read ring sugar");
            }
        }

    }


    public void getSugarPatterns(){
        getSugarChains();

        //patternToRing = new Hashtable<>();

        patternListLinearSugars = new ArrayList<>();
        for(IAtomContainer sugarAC : linearSugars){
            patternListLinearSugars.add(DfPattern.findSubstructure(sugarAC));
        }
        


    public IAtomContainer removeSugars(IAtomContainer molecule){
        //SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique );
        this.getSugarPatterns();

        IAtomContainer newMolecule = null;
        try {
            newMolecule = molecule.clone();

            //removing ring sugars

            int[][] g = GraphUtil.toAdjList(newMolecule);

            // efficient computation/partitioning of the ring systems
            RingSearch rs = new RingSearch(newMolecule, g);

            // isolated cycles don't need to be run
            List<IAtomContainer> isolatedRings = rs.isolatedRingFragments();
            for(IAtomContainer referenceRing : ringSugars){
                for(IAtomContainer isolatedRing : isolatedRings){
                    if (universalIsomorphismTester.isIsomorph(referenceRing, isolatedRing)) {
                        newMolecule.setProperty("CONTAINS_RING_SUGAR", 1);

                        //remove found ring
                        for(IAtom a : isolatedRing.atoms()){
                            newMolecule.removeAtom(a);
                        }
                    }
                }
            }


            //removing linear sugars
            for(DfPattern pattern : patternListLinearSugars){
                if (pattern.matches(molecule)) {

                    newMolecule.setProperty("CONTAINS_LINEAR_SUGAR", 1);
                    //remove sugar
                    // Mappings mappings = pattern.matchAll(newMolecule);

                    int[][] mappings = pattern.matchAll(newMolecule).uniqueAtoms().toArray();
                    //Mappings mappings = pattern.matchAll(newMolecule);
                    //Mappings mappings = pattern.matchAll(molecule);

                    for (int[] p : mappings) {

                        for (int i = 0; i < p.length; i++) {
                            IAtom atomToRemove = newMolecule.getAtom(i);
                            newMolecule.removeAtom(atomToRemove);
                        }
                    }

                }
            }
            //select only the biggest part of the molecule
            newMolecule = getBiggestComponent(newMolecule);

        } catch (CloneNotSupportedException | ConcurrentModificationException | IndexOutOfBoundsException  e) {
            e.printStackTrace();
            return null;
        } catch (CDKException e) {
            e.printStackTrace();
            return null;
        }

        return newMolecule;
    }


    IAtomContainer getBiggestComponent(IAtomContainer molecule){

        Map properties = molecule.getProperties();
        mcc = BeanUtil.getBean(MoleculeConnectivityChecker.class);
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
            if(nbheavyatoms<5){
                return null;
            }
        }
        else{

            return null;
        }
        molecule.setProperties(properties);
        return molecule;
    }






}
