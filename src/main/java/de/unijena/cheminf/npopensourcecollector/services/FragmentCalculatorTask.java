package de.unijena.cheminf.npopensourcecollector.services;


import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.mongocollections.Fragment;
import de.unijena.cheminf.npopensourcecollector.mongocollections.FragmentRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.annotation.Transient;
import org.springframework.stereotype.Component;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Propagation;
import org.springframework.transaction.annotation.Transactional;

import java.text.SimpleDateFormat;
import java.util.*;

@Component
@Transactional(propagation = Propagation.REQUIRED, readOnly = false)
public class FragmentCalculatorTask implements Runnable {

    @Autowired
    @Transient
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    @Transient
    FragmentRepository fragmentRepository;

    @Autowired
    @Transient
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;

    @Autowired
    @Transient
    SugarRemovalService sugarRemovalService;


    ElectronDonation model = ElectronDonation.cdk();
    CycleFinder cycles = Cycles.cdkAromaticSet();
    Aromaticity aromaticity = new Aromaticity(model, cycles);
    UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();




    private final int height = 2;

    List<UniqueNaturalProduct> batchOfNaturalProducts;

    Integer taskid;


    @Override
    public void run() {

        this.uniqueNaturalProductRepository = BeanUtil.getBean(UniqueNaturalProductRepository.class);
        this.fragmentRepository = BeanUtil.getBean(FragmentRepository.class);
        this.atomContainerToUniqueNaturalProductService = BeanUtil.getBean(AtomContainerToUniqueNaturalProductService.class);
        this.sugarRemovalService = BeanUtil.getBean(SugarRemovalService.class);
        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique); //Unique - canonical SMILES string, different atom ordering produces the same* SMILES. No isotope or stereochemistry encoded.


        System.out.println("Computing NP fragments for task "+taskid);
        for(UniqueNaturalProduct np : batchOfNaturalProducts){

            Double npl_score = 0.0;
            Double npl_score_with_sugar = 0.0;
            Double npl_score_noh = 0.0;


            IAtomContainer acFull = atomContainerToUniqueNaturalProductService.createAtomContainer(np);

            IAtomContainer acSugarFree = this.sugarRemovalService.removeSugars(acFull);

            if(acSugarFree != null) {

                //molecule contains sugar and is not only sugar OR does not contain sugar
                try {
                    if(!universalIsomorphismTester.isIsomorph(acSugarFree, acFull)) {
                        np.setContains_sugar(1);
                        np.setContains_ring_sugars(acSugarFree.getProperty("CONTAINS_RING_SUGAR"));
                        np.setContains_linear_sugars(acSugarFree.getProperty("CONTAINS_LINEAR_SUGAR"));
                    }
                    else{
                        //molecule doesn't contain sugar
                        np.setContains_sugar(0);
                    }
                } catch (CDKException e) {
                    System.out.println("Failed detecting isomorphism between molecules with and without sugar");
                    e.printStackTrace();
                }

                //counting atoms for sugar free molecule
                np.setSugar_free_total_atom_number(acSugarFree.getAtomCount());
                np.setSugar_free_heavy_atom_number(computeNumberOfHeavyAtoms(acSugarFree));


                //fragmenting the 2 versions of the molecule
                Hashtable<String, Integer> fragmentsWithSugar = generateCountedAtomSignatures(acFull, height);
                Hashtable<String, Integer> fragmentsWithoutSugar = generateCountedAtomSignatures(acSugarFree, height);



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
                        newFragment.presentInMolecules.add(np);
                        foundFragment = fragmentRepository.save(newFragment);
                    }
                    else{
                        foundFragment.presentInMolecules.add(np);
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
                        newFragment.presentInMolecules.add(np);
                        foundFragment = fragmentRepository.save(newFragment);
                    }
                    else{
                        foundFragment.presentInMolecules.add(np);
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


                try {
                    IAtomContainer nm = AtomContainerManipulator.removeHydrogens(acSugarFree);
                    np.setSugar_free_smiles(smilesGenerator.create(nm));
                } catch (CDKException e) {
                    e.printStackTrace();
                }


                uniqueNaturalProductRepository.save(np);
                //System.out.println("Saved natural products without sugar and NPL score");

            }
            else{
                //molecule is only sugar - need to work only on the sugar version

                //molecule after sugar removal is null - only sugar
                np.setContains_sugar(2);

                Hashtable<String, Integer> fragmentsWithSugar = generateCountedAtomSignatures(acFull, height);


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
                        newFragment.presentInMolecules.add(np);
                        foundFragment = fragmentRepository.save(newFragment);
                    }

                    npl_score_with_sugar = npl_score_with_sugar + (foundFragment.getScorenp() * fragmentsWithSugar.get(f));

                    np.addFragmentWithSugar(f, fragmentsWithSugar.get(f));

                }
                npl_score_with_sugar = npl_score_with_sugar / np.getTotal_atom_number();
                np.setNpl_sugar_score(npl_score_with_sugar);


                uniqueNaturalProductRepository.save(np);

                //System.out.println("Saved natural products with only sugars and NPL score");
            }


            //Calculate Murko Framework

            try {
                MurckoFragmenter murckoFragmenter = new MurckoFragmenter(true, 3);
                murckoFragmenter.generateFragments(acFull);
                if(murckoFragmenter.getFragments() != null && murckoFragmenter.getFragments().length >0) {
                    np.setMurko_framework(murckoFragmenter.getFrameworks()[0]);
                    uniqueNaturalProductRepository.save(np);
                }

            } catch (CDKException | NullPointerException e) {
                //e.printStackTrace();
                System.out.println("Failed creating Murcko fragment");
                np.setMurko_framework("");
                uniqueNaturalProductRepository.save(np);
            }


        }
        Date date = new Date();
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
        System.out.println("Task "+taskid+" finished at: "+formatter.format(date)+"\n");


    }

    public List<UniqueNaturalProduct> getBatchOfNaturalProducts() {
        return batchOfNaturalProducts;
    }

    public void setBatchOfNaturalProducts(List<UniqueNaturalProduct> batchOfNaturalProducts) {
        this.batchOfNaturalProducts = batchOfNaturalProducts;
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




}
