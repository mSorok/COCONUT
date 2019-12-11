package de.unijena.cheminf.npopensourcecollector.services;


import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.misc.MoleculeConnectivityChecker;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.BondManipulator;
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

        linearSugars = new ArrayList<>();
        ringSugars = new ArrayList<>();
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String[] linearSugarsSmilesList = {"C(C(C(C(C(C=O)O)O)O)O)O", "C(C(CC(C(CO)O)O)O)(O)=O", "C(C(C(CC(=O)O)O)O)O",
                "C(C(C(C(C(CO)O)O)O)=O)O", "C(C(C(C(C(CO)O)O)O)O)O", "C(C(C(C(CC=O)O)O)O)O",
                "OCC(O)C(O)C(O)C(O)CO", "O=CC(O)C(O)C(O)C(O)CO",  "CCCCC(O)C(=O)O", "CC(=O)CC(=O)CCC(=O)O",
                "O=C(O)CC(O)CC(=O)O", "O=C(O)C(=O)C(=O)C(O)C(O)CO",
                "O=C(O)CCC(O)C(=O)O", "O=CC(O)C(O)C(O)C(O)CO", "O=C(CO)C(O)C(O)CO"};

        String [] ringSugarsSmilesList = { "C1CCOC1", "C1CCOCC1", "C1CCCOCC1"};



        //adding linear sugars to list
        for (String smiles : linearSugarsSmilesList) {
            try {
                linearSugars.add(smilesParser.parseSmiles(smiles));
            } catch (InvalidSmilesException ex) {
                System.out.println("could not read linear sugar");
            }
        }


        //adding ring sugars to list
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

        patternListLinearSugars = new ArrayList<>();
        for(IAtomContainer sugarAC : linearSugars){
            patternListLinearSugars.add(DfPattern.findSubstructure(sugarAC));
        }

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
            List<IAtomContainer> isolatedRings = rs.isolatedRingFragments();//
            List<IAtomContainer> fusedRings = rs.fusedRingFragments();//
            for(IAtomContainer referenceRing : ringSugars){//
                for(IAtomContainer isolatedRing : isolatedRings){//
                    if (universalIsomorphismTester.isIsomorph(referenceRing, isolatedRing)) {
                        boolean tmpAreAllExocyclicBondsSingle = areAllExocyclicBondsSingle(isolatedRing, newMolecule);
                        if (!tmpAreAllExocyclicBondsSingle) {
                            continue;
                        }
                        int tmpExocyclicOxygenCount = getAttachedOxygenAtomCount(isolatedRing, newMolecule);
                        // To do: Count only heavy atoms? But normally, Hs should not be included in the atom containers using RingSearch
                        int tmpAtomsInRing = isolatedRing.getAtomCount();
                        boolean tmpAreEnoughOxygensAttached = doesRingHaveEnoughOxygenAtomsAttached(tmpAtomsInRing, tmpExocyclicOxygenCount);
                        if (!tmpAreEnoughOxygensAttached) {
                            continue;
                        }
                        //</Jonas>
                        newMolecule.setProperty("CONTAINS_RING_SUGAR", 1);

                        //remove found ring
                        for(IAtom a : isolatedRing.atoms()){
                            //<Jonas>
                            if (!newMolecule.contains(a)) {
                                continue;
                            }
                            //</Jonas>
                            newMolecule.removeAtom(a);
                        }
                    }
                }
            }


            //removing linear sugars
            boolean containsLinearSugar = false;
            for(DfPattern pattern : patternListLinearSugars){
                //if (pattern.matches(molecule)) {
                //boolean patternPartOfARing = false;

                //remove sugar

                int[][] mappings = pattern.matchAll(newMolecule).uniqueAtoms().toArray();

                for (int[] p : mappings) {
                    boolean matchPartOfRing = false;
                    for (int i = 0; i < p.length; i++) {
                        if( newMolecule.getAtomCount()>= i ){
                            for(IAtomContainer fr : fusedRings){
                                if(fr.contains(newMolecule.getAtom(i))){
                                    //mapped atom is in a ring
                                    //patternPartOfARing = true;
                                    matchPartOfRing = true;
                                    break;
                                }
                            }
                            if (matchPartOfRing) {
                                break;
                            }
                            for(IAtomContainer ir : isolatedRings){
                                if(ir.contains(newMolecule.getAtom(i))){
                                    //mapped atom is in a ring
                                    //patternPartOfARing = true;
                                    matchPartOfRing = true;
                                    break;
                                }
                            }
                            if (matchPartOfRing) {
                                break;
                            }

                        }

                    }
                    if(!matchPartOfRing ){
                        //for (int[] p : mappings) {
                        ArrayList<IAtom> atomsToRemove = new ArrayList<>();

                        for (int i = 0; i < p.length; i++) {
                            if( newMolecule.getAtomCount()>= i ) {
                                //<Jonas> added try catch block
                                try {
                                    IAtom atomToRemove = newMolecule.getAtom(i);
                                    atomsToRemove.add(atomToRemove);
                                } catch (IndexOutOfBoundsException anException) {
                                    anException.printStackTrace();
                                }
                                //</Jonas>
                            }
                        }
                        if(atomsToRemove.size()>=4){
                            containsLinearSugar = true;
                        }
                        for(IAtom a : atomsToRemove){
                            //<Jonas>
                            if (!newMolecule.contains(a)) {
                                continue;
                            }
                            //</Jonas>
                            newMolecule.removeAtom(a);
                        }

                        //}

                    }
                }


            }
            if(containsLinearSugar){
                newMolecule.setProperty("CONTAINS_LINEAR_SUGAR", 1);
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
            if(nbheavyatoms<=6){
                return null;
            }
        }
        else{

            return null;
        }
        molecule.setProperties(properties);
        return molecule;
    }





    /**
     * Checks whether all exocyclic bonds connected to a given ring fragment of an original atom container are of single
     * order.
     * <p>
     *     The method iterates over all cyclic atoms and all of their bonds. So the runtime scales linear with the number
     *     of cyclic atoms and their connected bonds.
     * </p>
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return true, if all exocyclic bonds connected to the ring are of single order
     * @throws NullPointerException if a parameter is 'null'
     * @author Jonas
     */
    public static boolean areAllExocyclicBondsSingle(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule) throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule is 'null'");
        int tmpAtomCountInRing = aRingToTest.getAtomCount();
        int tmpArrayListInitCapacity = tmpAtomCountInRing * 3;
        List<IBond> tmpExocyclicBondsList = new ArrayList<>(tmpArrayListInitCapacity);
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            // To do: 'NoSuchAtomException: Atom is not a member of this AtomContainer' is thrown sometimes in the following line
            List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpRingAtom);
            for (IBond tmpBond : tmpConnectedBondsList) {
                boolean tmpIsInRing = aRingToTest.contains(tmpBond);
                if (!tmpIsInRing) {
                    tmpExocyclicBondsList.add(tmpBond);
                }
            }
        }
        return (BondManipulator.getMaximumBondOrder(tmpExocyclicBondsList) == IBond.Order.SINGLE);
    }

    // To do: Also count N and S (or all hetero atoms)?
    /**
     * Gives the number of attached exocyclic oxygen atoms of a given ring fragment of an original atom container.
     * <p>
     *     The method iterates over all cyclic atoms and all of their connected atoms. So the runtime scales linear with the number
     *     of cyclic atoms and their connected atoms.
     * </p>
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return number of attached exocyclic oxygen atoms
     * @throws NullPointerException if a parameter is 'null'
     * @author Jonas
     */
    public static int getAttachedOxygenAtomCount(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule) throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule is 'null'");
        int tmpExocyclicOxygenCounter = 0;
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            // To do: 'NoSuchAtomException: Atom is not a member of this AtomContainer' is thrown sometimes in the following line
            List<IAtom> tmpConnectedAtomsList = anOriginalMolecule.getConnectedAtomsList(tmpRingAtom);
            for (IAtom tmpConnectedAtom : tmpConnectedAtomsList) {
                String tmpSymbol = tmpConnectedAtom.getSymbol();
                boolean tmpIsOxygen = tmpSymbol.matches("O");
                boolean tmpIsInRing = aRingToTest.contains(tmpConnectedAtom);
                if (tmpIsOxygen && !tmpIsInRing) {
                    tmpExocyclicOxygenCounter++;
                }
            }
        }
        return tmpExocyclicOxygenCounter;
    }

    /**
     * Simple decision making function for deciding whether a possible sugar ring has enough exocyclic oxygen atom
     * attached to it. This number should be higher or equal to the number of atoms in the ring (including the cyclic
     * oxygen atom) divided by two. So at least 3 attached exocyclic oxygen atoms for a six-membered ring, 2 for a
     * five-membered ring etc. A simple integer division is used.
     *
     * @param aNumberOfAtomsInRing number of atoms in the possible sugar ring, including the cyclic oxygen atom
     * @param aNumberOfAttachedExocyclicOxygenAtoms number of attached exocyclic oxygen atom of the ring under
     *                                              investigation
     * @return true, if the number of attached exocyclic oxygen atoms is at least half of the number of atoms in the ring
     */
    public static boolean doesRingHaveEnoughOxygenAtomsAttached(int aNumberOfAtomsInRing, int aNumberOfAttachedExocyclicOxygenAtoms) {
        return (aNumberOfAttachedExocyclicOxygenAtoms >= (aNumberOfAtomsInRing / 2));
    }




}
