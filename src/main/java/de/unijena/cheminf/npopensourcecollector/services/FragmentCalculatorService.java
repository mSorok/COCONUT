package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.misc.LinearSugars;
import de.unijena.cheminf.npopensourcecollector.mongocollections.Fragment;
import de.unijena.cheminf.npopensourcecollector.mongocollections.FragmentRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

@Service
public class FragmentCalculatorService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    FragmentRepository fragmentRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;


    ElectronDonation model = ElectronDonation.cdk();
    CycleFinder cycles = Cycles.cdkAromaticSet();
    Aromaticity aromaticity = new Aromaticity(model, cycles);


    private final LinearSugars linearSugarChains = LinearSugars.getInstance();

    private final int height = 2;



    public void doWork(){

        System.out.println("Start fragmenting natural products");


        List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();

        int count=1;
        int total=allNP.size();

        for(UniqueNaturalProduct np : allNP){


            IAtomContainer acFull = atomContainerToUniqueNaturalProductService.createAtomContainer(np);

            IAtomContainer acSugarFree = removeSugars(acFull);

            if(acSugarFree != null && acSugarFree.getAtomCount()>0) {

                //counting atoms for sugar free molecule
                np.setSugar_free_total_atom_number(acSugarFree.getAtomCount());
                np.setSugar_free_heavy_atom_number(computeNumberOfHeavyAtoms(acSugarFree));


                //fragmenting the 2 versions of the molecule
                Hashtable<String, Integer> fragmentsWithSugar = generateCountedAtomSignatures(acFull, height);
                Hashtable<String, Integer> fragmentsWithoutSugar = generateCountedAtomSignatures(acSugarFree, height);

                Double npl_score = 0.0;
                Double npl_score_with_sugar = 0.0;
                Double npl_score_noh = 0.0;



                //computing the NPL score with the Sugar
                for (String f : fragmentsWithSugar.keySet()) {

                    Fragment foundFragment = fragmentRepository.findBySignatureAndWithsugar(f, 1);

                    if(foundFragment==null){
                        //it is a new fragment!
                        Fragment newFragment = new Fragment();
                        newFragment.setHeight(height);
                        newFragment.setWith_sugar(1);
                        newFragment.setSignature(f);
                        newFragment.setScorenp(1.0);
                        foundFragment = fragmentRepository.save(newFragment);
                    }

                    npl_score_with_sugar = npl_score_with_sugar + (foundFragment.getScorenp() * fragmentsWithSugar.get(f));

                    np.addFragmentWithSugar(f, fragmentsWithSugar.get(f));

                }
                npl_score_with_sugar = npl_score_with_sugar / np.getTotal_atom_number();
                np.setNpl_sugar_score(npl_score_with_sugar);

                //Computing the NPL score without the sugar
                for (String f : fragmentsWithoutSugar.keySet()) {
                    Fragment foundFragment = fragmentRepository.findBySignatureAndWithsugar(f, 0);

                    if(foundFragment==null){
                        //it is a new fragment!
                        Fragment newFragment = new Fragment();
                        newFragment.setHeight(height);
                        newFragment.setWith_sugar(0);
                        newFragment.setSignature(f);
                        newFragment.setScorenp(1.0);
                        foundFragment = fragmentRepository.save(newFragment);
                    }

                    npl_score = npl_score + (foundFragment.getScorenp() * fragmentsWithoutSugar.get(f));

                    np.addFragment(f, fragmentsWithoutSugar.get(f) );

                    //For the score without fragments starting by a H
                    if (!f.startsWith("[H]")) {
                        npl_score_noh = npl_score_noh + (foundFragment.getScorenp() * fragmentsWithoutSugar.get(f));
                    }
                }

                npl_score = npl_score / np.getSugar_free_total_atom_number();
                np.setNpl_score(npl_score);

                npl_score_noh = npl_score_noh / np.getSugar_free_heavy_atom_number();
                np.setNpl_noh_score(npl_score_noh);

                uniqueNaturalProductRepository.save(np);
            }
            count++;
            if(count%10000==0){
                System.out.println("Molecules fragmented: "+count+" ("+(double)count/(double)total+"% )");

            }
        }


        System.out.println("Done fragmenting natural products");
    }





    public List<String> generateAtomSignatures(IAtomContainer atomContainer, Integer height) {

        List<String> atomSignatures = new ArrayList<>();



        //atomContainer = calculateAromaticity(atomContainer);

        if( atomContainer != null  && !atomContainer.isEmpty()) {

            for (IAtom atom : atomContainer.atoms()) {
                try {
                    AtomSignature atomSignature = new AtomSignature(atom, height, atomContainer);
                    atomSignatures.add(atomSignature.toCanonicalString());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            return atomSignatures;
        }
        else{
            return null;
        }
    }


    public Hashtable<String, Integer> generateCountedAtomSignatures(IAtomContainer atomContainer, Integer height) {

        List<String> atomSignatures = new ArrayList<>();

        Hashtable<String, Integer> countedAtomSignatures = new Hashtable<>();



        //atomContainer = calculateAromaticity(atomContainer);

        if(atomContainer !=null && !atomContainer.isEmpty()) {

            for (IAtom atom : atomContainer.atoms()) {
                try {
                    AtomSignature atomSignature = new AtomSignature(atom, height, atomContainer);
                    atomSignatures.add(atomSignature.toCanonicalString());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            for(String signature : atomSignatures){
                if(countedAtomSignatures.containsKey(signature)){
                    countedAtomSignatures.put(signature, countedAtomSignatures.get(signature)+1);
                }
                else{
                    countedAtomSignatures.put(signature,1);
                }

            }


            return countedAtomSignatures;
        }
        else{
            return null;
        }
    }





    private IAtomContainer removeSugars(IAtomContainer molecule){

        IAtomContainer newMolecule = null;
        try {
            newMolecule = molecule.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }

        try {

            IRingSet ringset = Cycles.sssr(newMolecule).toRingSet();

            // RING SUGARS
            for (IAtomContainer one_ring : ringset.atomContainers()) {
                try {
                    IMolecularFormula molecularFormula = MolecularFormulaManipulator.getMolecularFormula(one_ring);
                    String formula = MolecularFormulaManipulator.getString(molecularFormula);
                    IBond.Order bondorder = AtomContainerManipulator.getMaximumBondOrder(one_ring);

                    if (formula.equals("C5O") | formula.equals("C4O") | formula.equals("C6O")) {
                        if (IBond.Order.SINGLE.equals(bondorder)) {
                            if (shouldRemoveRing(one_ring, newMolecule, ringset) == true) {
                                for (IAtom atom : one_ring.atoms()) {
                                    {

                                        newMolecule.removeAtom(atom);
                                    }
                                }
                            }

                        }
                    }
                }catch(NullPointerException e){
                    return null;
                }
            }
            Map<Object, Object> properties = newMolecule.getProperties();
            IAtomContainerSet molset = ConnectivityChecker.partitionIntoMolecules(newMolecule);
            for (int i = 0; i < molset.getAtomContainerCount(); i++) {
                molset.getAtomContainer(i).setProperties(properties);
                int size = molset.getAtomContainer(i).getBondCount();
                if (size >= 5) {
                    if (!linearSugarChains.hasSugarChains(molset.getAtomContainer(i), ringset.getAtomContainerCount())) {

                        return (IAtomContainer) molset.getAtomContainer(i);
                    }
                }
            }
            //
        } catch (NullPointerException e) {
        } catch (CDKException e) {
        }
        return null;

    }





    private boolean shouldRemoveRing(IAtomContainer possibleSugarRing, IAtomContainer molecule, IRingSet sugarRingsSet) {

        boolean shouldRemoveRing = false;
        List<IAtom> allConnectedAtoms = new ArrayList<IAtom>();
        List<IBond> bonds = new ArrayList<IBond>();
        int oxygenAtomCount = 0;

        IRingSet connectedRings = sugarRingsSet.getConnectedRings((IRing) possibleSugarRing);

        /*
         * get bonds to check for bond order of connected atoms in a sugar ring
         *
         */
        for (IAtom atom : possibleSugarRing.atoms()) {
            bonds.addAll(molecule.getConnectedBondsList(atom));
        }

        if (IBond.Order.SINGLE.equals(BondManipulator.getMaximumBondOrder(bonds))
                && connectedRings.getAtomContainerCount() == 0) {

            /*
             * get connected atoms of all atoms in sugar ring to check for glycoside bond
             */
            for (IAtom atom : possibleSugarRing.atoms()) {
                List<IAtom> connectedAtoms = molecule.getConnectedAtomsList(atom);
                allConnectedAtoms.addAll(connectedAtoms);
            }

            for (IAtom connected_atom : allConnectedAtoms) {
                if (!possibleSugarRing.contains(connected_atom)) {
                    if (connected_atom.getSymbol().matches("O")) {
                        oxygenAtomCount++;
                    }
                }
            }
            if (oxygenAtomCount > 0) {
                return true;
            }
        }
        return shouldRemoveRing;
    }


    public Integer computeNumberOfHeavyAtoms(IAtomContainer ac){
        Integer numberHeavyAtoms = 0;
        for(IAtom a : ac.atoms() ){
            if(!a.getSymbol().equals("H")){
                numberHeavyAtoms++;
            }
        }
        return numberHeavyAtoms;
    }



}
