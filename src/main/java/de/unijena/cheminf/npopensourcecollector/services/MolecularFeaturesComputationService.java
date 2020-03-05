package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SyntheticMolecule;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SyntheticMoleculeRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.*;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinderUtility;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;



@Service
public class MolecularFeaturesComputationService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;

    @Autowired
    SyntheticMoleculeRepository syntheticMoleculeRepository;

    @Autowired
    AtomContainerToSyntheticMoleculeService atomContainerToSyntheticMoleculeService;

    PubchemFingerprinter pubchemFingerprinter = new PubchemFingerprinter( SilentChemObjectBuilder.getInstance() );

    CircularFingerprinter circularFingerprinter = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);

    KlekotaRothFingerprinter klekotaRothFingerprinter = new KlekotaRothFingerprinter();

    HybridizationFingerprinter hybridizationFingerprinter = new HybridizationFingerprinter();

    MACCSFingerprinter maccsFingerprinter = new MACCSFingerprinter();

    ShortestPathFingerprinter shortestPathFingerprinter = new ShortestPathFingerprinter();

    SubstructureFingerprinter substructureFingerprinter = new SubstructureFingerprinter();

    Aromaticity aromaticityModel = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
    ErtlFunctionalGroupsFinder ertlFunctionalGroupsFinder  = ErtlFunctionalGroupsFinderUtility.getErtlFunctionalGroupsFinderGeneralizingMode();
    MoleculeHashGenerator efgHashGenerator = ErtlFunctionalGroupsFinderUtility.getFunctionalGroupHashGenerator();
    SmilesGenerator efgSmilesGenerator = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);






    public void doWorkRecompute(){
        System.out.println("Calculating additional features for unique molecules (only incomplete)");

        List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAllByApolComputed();
        for(UniqueNaturalProduct np : allNP){

            np = computeFeatures(np);
            np = computeFingerprints(np);

            uniqueNaturalProductRepository.save(np);
        }
        System.out.println("done");
    }



    public void doWork(){
        System.out.println("Calculating additional features for unique molecules");


        List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();

        for(UniqueNaturalProduct np : allNP){

            np = computeFeatures(np);
            np = computeFingerprints(np);
            np = computeErtlFunctionalGroups(np);

            uniqueNaturalProductRepository.save(np);
        }



        System.out.println("done");
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

            np.setPubchemFingerprint(pubchemFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setCircularFingerprint(circularFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setKlekotaRothFingerprint(klekotaRothFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setHybridizationFingerprint(hybridizationFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setMaccsFingerprint(maccsFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setShortestPathFingerprint(shortestPathFingerprinter.getBitFingerprint(ac).asBitSet());
            np.setSubstructureFingerprint(substructureFingerprinter.getBitFingerprint(ac).asBitSet());

        } catch (CDKException | UnsupportedOperationException e) {
            e.printStackTrace();
        }
        return np;
    }


    public UniqueNaturalProduct computeErtlFunctionalGroups(UniqueNaturalProduct np){

        IAtomContainer ac = atomContainerToUniqueNaturalProductService.createAtomContainer(np);
        List<IAtomContainer> functionalGroupsGeneralized;

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


        ac = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(ac, aromaticityModel);
        functionalGroupsGeneralized = ertlFunctionalGroupsFinder.find(ac, false);
        if (!functionalGroupsGeneralized.isEmpty()) {

            np.ertlFunctionalFragments = new Hashtable<>();
            np.ertlFunctionalFragmentsPseudoSmiles = new Hashtable<>();

            HashMap<Long, IAtomContainer> tmpResultsMap = new HashMap<>(functionalGroupsGeneralized.size(), 1);
            for (IAtomContainer functionalGroup : functionalGroupsGeneralized) {
                Long hashCode = efgHashGenerator.generate(functionalGroup);
                if (tmpResultsMap.keySet().contains(hashCode)) {
                    int tmpFrequency = tmpResultsMap.get(hashCode).getProperty("FREQUENCY");
                    tmpResultsMap.get(hashCode).setProperty("FREQUENCY", tmpFrequency + 1);
                } else {
                    functionalGroup.setProperty("FREQUENCY", 1);
                    tmpResultsMap.put(hashCode, functionalGroup);
                }
            }


            for (Long tmpHashCode : tmpResultsMap.keySet()) {
                IAtomContainer tmpFunctionalGroup = tmpResultsMap.get(tmpHashCode);
                String tmpFGSmilesCode = null;
                try {
                    tmpFGSmilesCode = efgSmilesGenerator.create(tmpFunctionalGroup);

                    String tmpFGPseudoSmilesCode = ErtlFunctionalGroupsFinderUtility.createPseudoSmilesCode(tmpFunctionalGroup);


                    int tmpFrequency = tmpFunctionalGroup.getProperty("FREQUENCY");

                    np.ertlFunctionalFragments.put(tmpFGSmilesCode, tmpFrequency);
                    np.ertlFunctionalFragmentsPseudoSmiles.put(tmpFGPseudoSmilesCode, tmpFrequency);

                } catch (CDKException e) {
                    e.printStackTrace();
                }
            }

        }




        return np;

    }


    public void predictStereochemistry(){
        //TODO - eventually one day, but very high combinatorics
    }



    public void doWorkForSM(){

        System.out.println("Calculating additional features for synthetic molecules");

        List<SyntheticMolecule> allSM = syntheticMoleculeRepository.findAll();

        //Constructors for descriptors
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        AtomContainerManipulator acm = new AtomContainerManipulator();



        for(SyntheticMolecule sm : allSM){
            IAtomContainer acFull = atomContainerToSyntheticMoleculeService.createAtomContainer(sm);

            //get bond count
            IBond[] bonds = acm.getBondArray(acFull);
            int bondCount = 0;
            for(IBond b : bonds){
                if(b.getAtomCount() == 2){
                    if(!b.getAtom(0).getSymbol().equals("H") && !b.getAtom(1).getSymbol().equals("H")){
                        bondCount++;
                    }
                }
            }
            sm.setBond_count(bondCount);

            try{

                //AlogP
                ALOGPDescriptor alogpDescriptor = new ALOGPDescriptor();
                alogpDescriptor.initialise(builder);
                DescriptorValue alogpvalue = alogpDescriptor.calculate(acFull);
                DoubleArrayResult alogpresults = (DoubleArrayResult) alogpvalue.getValue();
                sm.setAlogp(alogpresults.get(0));
                sm.setAlogp2(alogpresults.get(1));
                sm.setAmralogp(alogpresults.get(2));


                APolDescriptor aPolDescriptor = new APolDescriptor();
                aPolDescriptor.initialise(builder);
                DescriptorValue apolvalue = aPolDescriptor.calculate(acFull);
                DoubleResult apolresult = (DoubleResult) apolvalue.getValue();
                sm.setApol(apolresult.doubleValue());

                try {
                    BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
                    bcutDescriptor.initialise(builder);
                    sm.bcutDescriptor = new ArrayList<>();
                    DescriptorValue bcutvalue = bcutDescriptor.calculate(acFull);
                    DoubleArrayResult bcutResults = (DoubleArrayResult) bcutvalue.getValue();
                    for (int i = 0; i < 6; i++) {
                        sm.bcutDescriptor.add(bcutResults.get(i));
                    }
                }catch(ArrayIndexOutOfBoundsException e){
                    e.printStackTrace();
                }

                BPolDescriptor bPolDescriptor = new BPolDescriptor();
                bPolDescriptor.initialise(builder);
                DescriptorValue bpolvalue = bPolDescriptor.calculate(acFull);
                DoubleResult bpolresult = (DoubleResult) bpolvalue.getValue();
                sm.setBpol(bpolresult.doubleValue());

                EccentricConnectivityIndexDescriptor eccentricConnectivityIndexDescriptor = new EccentricConnectivityIndexDescriptor();
                eccentricConnectivityIndexDescriptor.initialise(builder);
                DescriptorValue eccenValue = eccentricConnectivityIndexDescriptor.calculate(acFull);
                IntegerResult eccenResult = (IntegerResult) eccenValue.getValue();
                sm.setEccentricConnectivityIndexDescriptor(eccenResult.intValue());

                FMFDescriptor fmfDescriptor = new FMFDescriptor();
                fmfDescriptor.initialise(builder);
                DescriptorValue fmfValue = fmfDescriptor.calculate(acFull);
                DoubleResult fmfResult = (DoubleResult) fmfValue.getValue();
                sm.setFmfDescriptor(fmfResult.doubleValue());

                FractionalCSP3Descriptor fractionalCSP3Descriptor = new FractionalCSP3Descriptor();
                fractionalCSP3Descriptor.initialise(builder);
                DescriptorValue fsp3value = fractionalCSP3Descriptor.calculate(acFull);
                DoubleResult fsp3result = (DoubleResult) fsp3value.getValue();
                sm.setFsp3(fsp3result.doubleValue());

                FragmentComplexityDescriptor fragmentComplexityDescriptor = new FragmentComplexityDescriptor();
                fragmentComplexityDescriptor.initialise(builder);
                DescriptorValue fcdValue = fragmentComplexityDescriptor.calculate(acFull);
                DoubleResult fcdResults = (DoubleResult) fcdValue.getValue();
                sm.setFragmentComplexityDescriptor(fcdResults.doubleValue());


                GravitationalIndexDescriptor gravitationalIndexDescriptor = new GravitationalIndexDescriptor();
                gravitationalIndexDescriptor.initialise(builder);
                DescriptorValue gravValue = gravitationalIndexDescriptor.calculate(acFull);
                DoubleArrayResult gravResults = (DoubleArrayResult) gravValue.getValue();
                sm.setGravitationalIndexHeavyAtoms(gravResults.get(0));

                HBondAcceptorCountDescriptor hBondAcceptorCountDescriptor = new HBondAcceptorCountDescriptor();
                hBondAcceptorCountDescriptor.initialise(builder);
                DescriptorValue nhbaccValue = hBondAcceptorCountDescriptor.calculate(acFull);
                IntegerResult nhbaccResult = (IntegerResult) nhbaccValue.getValue();
                sm.sethBondAcceptorCount(nhbaccResult.intValue());

                HBondDonorCountDescriptor hBondDonorCountDescriptor = new HBondDonorCountDescriptor();
                hBondDonorCountDescriptor.initialise(builder);
                DescriptorValue nhbdonValue = hBondDonorCountDescriptor.calculate(acFull);
                IntegerResult nhbdonResult = (IntegerResult) nhbdonValue.getValue();
                sm.sethBondDonorCount(nhbdonResult.intValue());

                HybridizationRatioDescriptor hybridizationRatioDescriptor = new HybridizationRatioDescriptor();
                hybridizationRatioDescriptor.initialise(builder);
                DescriptorValue hybridRatioValue = hybridizationRatioDescriptor.calculate(acFull);
                DoubleResult hybridRationResult = (DoubleResult) hybridRatioValue.getValue();
                sm.setHybridizationRatioDescriptor(hybridRationResult.doubleValue());

                KappaShapeIndicesDescriptor kappaShapeIndicesDescriptor = new KappaShapeIndicesDescriptor();
                kappaShapeIndicesDescriptor.initialise(builder);
                DescriptorValue kappaShapeValues = kappaShapeIndicesDescriptor.calculate(acFull);
                DoubleArrayResult kappaShapeResults = (DoubleArrayResult) kappaShapeValues.getValue();
                sm.setKappaShapeIndex1(kappaShapeResults.get(0));
                sm.setKappaShapeIndex2(kappaShapeResults.get(1));
                sm.setKappaShapeIndex3(kappaShapeResults.get(2));

                MannholdLogPDescriptor mannholdLogPDescriptor = new MannholdLogPDescriptor();
                mannholdLogPDescriptor.initialise(builder);
                DescriptorValue manholdLogpValues = mannholdLogPDescriptor.calculate(acFull);
                DoubleResult manholdLogpResult = (DoubleResult) manholdLogpValues.getValue();
                sm.setManholdlogp(manholdLogpResult.doubleValue());

                PetitjeanNumberDescriptor petitjeanNumberDescriptor = new PetitjeanNumberDescriptor();
                petitjeanNumberDescriptor.initialise(builder);
                DescriptorValue petitjeanNumnerValue = petitjeanNumberDescriptor.calculate(acFull);
                DoubleResult petitjeanResult = (DoubleResult) petitjeanNumnerValue.getValue();
                sm.setPetitjeanNumber(petitjeanResult.doubleValue());

                PetitjeanShapeIndexDescriptor petitjeanShapeIndexDescriptor = new PetitjeanShapeIndexDescriptor();
                petitjeanShapeIndexDescriptor.initialise(builder);
                DescriptorValue petitjeanShapeValues = petitjeanShapeIndexDescriptor.calculate(acFull);
                DoubleArrayResult petitjeanShapeResults = (DoubleArrayResult) petitjeanShapeValues.getValue();
                sm.setPetitjeanShapeTopo(petitjeanShapeResults.get(0));
                sm.setPetitjeanShapeGeom(petitjeanShapeResults.get(1));

                RuleOfFiveDescriptor ruleOfFiveDescriptor = new RuleOfFiveDescriptor();
                ruleOfFiveDescriptor.initialise(builder);
                DescriptorValue ruleOfFiveValue = ruleOfFiveDescriptor.calculate(acFull);
                IntegerResult ruleOfFiveResult = (IntegerResult) ruleOfFiveValue.getValue();
                sm.setLipinskiRuleOf5Failures(ruleOfFiveResult.intValue());

                /*SmallRingDescriptor smallRingDescriptor = new SmallRingDescriptor();
                smallRingDescriptor.initialise(builder);
                DescriptorValue smallRingValue = smallRingDescriptor.calculate(acFull);
                IntegerResult smallRingResult = (IntegerResult) smallRingValue.getValue();
                np.setNumberSmallRingsDescriptor(smallRingResult.intValue());*/


                SpiroAtomCountDescriptor spiroAtomCountDescriptor = new SpiroAtomCountDescriptor();
                spiroAtomCountDescriptor.initialise(builder);
                DescriptorValue spiroatomValue = spiroAtomCountDescriptor.calculate(acFull);
                IntegerResult spiroatomResult = (IntegerResult) spiroatomValue.getValue();
                sm.setNumberSpiroAtoms(spiroatomResult.intValue());

                VABCDescriptor vabcDescriptor = new VABCDescriptor();
                vabcDescriptor.initialise(builder);
                DescriptorValue vabcValue = vabcDescriptor.calculate(acFull);
                DoubleResult vabcResult = (DoubleResult) vabcValue.getValue();
                sm.setVabcDescriptor(vabcResult.doubleValue());

                VAdjMaDescriptor vAdjMaDescriptor = new VAdjMaDescriptor();
                vAdjMaDescriptor.initialise(builder);
                DescriptorValue vadjmaValue = vAdjMaDescriptor.calculate(acFull);
                DoubleResult vadjmaResult = (DoubleResult) vadjmaValue.getValue();
                sm.setVertexAdjMagnitude(vadjmaResult.doubleValue());

                WienerNumbersDescriptor wienerNumbersDescriptor = new WienerNumbersDescriptor();
                wienerNumbersDescriptor.initialise(builder);
                DescriptorValue wienerNumbersValue = wienerNumbersDescriptor.calculate(acFull);
                DoubleArrayResult wienerNumbersResult = (DoubleArrayResult) wienerNumbersValue.getValue();
                sm.setWeinerPathNumber(wienerNumbersResult.get(0));
                sm.setWeinerPolarityNumber(wienerNumbersResult.get(1));

                XLogPDescriptor xLogPDescriptor = new XLogPDescriptor();
                xLogPDescriptor.initialise(builder);
                DescriptorValue xlogpValues = xLogPDescriptor.calculate(acFull);
                DoubleResult xlogpResult = (DoubleResult) xlogpValues.getValue();
                sm.setXlogp(xlogpResult.doubleValue());

                ZagrebIndexDescriptor zagrebIndexDescriptor = new ZagrebIndexDescriptor();
                zagrebIndexDescriptor.initialise(builder);
                DescriptorValue zagrebIndexValue = zagrebIndexDescriptor.calculate(acFull);
                DoubleResult zagrebIndexresult = (DoubleResult) zagrebIndexValue.getValue();
                sm.setZagrebIndex(zagrebIndexresult.doubleValue());

                TPSADescriptor tpsaDescriptor = new TPSADescriptor();
                tpsaDescriptor.initialise(builder);
                DescriptorValue tpsaValue = tpsaDescriptor.calculate(acFull);
                DoubleResult tpsaResult = (DoubleResult) tpsaValue.getValue();
                sm.setTopoPSA(tpsaResult.doubleValue());

                FractionalPSADescriptor fractionalPSADescriptor = new FractionalPSADescriptor();
                fractionalPSADescriptor.initialise(builder);
                DescriptorValue ftpsaValue = fractionalPSADescriptor.calculate(acFull);
                DoubleResult ftpsaResult = (DoubleResult) ftpsaValue.getValue();
                sm.setTpsaEfficiency(ftpsaResult.doubleValue());


                //compute the clean smiles
                SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique); //Unique - canonical SMILES string, different atom ordering produces the same* SMILES. No isotope or stereochemistry encoded.
                try {
                    IAtomContainer nm = AtomContainerManipulator.removeHydrogens(acFull);
                    String cleanSmiles =  smilesGenerator.create(nm);
                    sm.setClean_smiles(cleanSmiles);
                } catch (CDKException e) {
                    e.printStackTrace();
                }





            } catch (CDKException | OutOfMemoryError e) {
                e.printStackTrace();
            }



            syntheticMoleculeRepository.save(sm);


        }


        System.out.println("done");

    }



    public UniqueNaturalProduct computeFeatures(UniqueNaturalProduct np){

        //Constructors for descriptors
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();

        IAtomContainer acFull = atomContainerToUniqueNaturalProductService.createAtomContainer(np);

        try {
            //AlogP
            ALOGPDescriptor alogpDescriptor = new ALOGPDescriptor();
            alogpDescriptor.initialise(builder);
            DescriptorValue alogpvalue = alogpDescriptor.calculate(acFull);
            DoubleArrayResult alogpresults = (DoubleArrayResult) alogpvalue.getValue();
            np.setAlogp(alogpresults.get(0));
            np.setAlogp2(alogpresults.get(1));
            np.setAmralogp(alogpresults.get(2));


            APolDescriptor aPolDescriptor = new APolDescriptor();
            aPolDescriptor.initialise(builder);
            DescriptorValue apolvalue = aPolDescriptor.calculate(acFull);
            DoubleResult apolresult = (DoubleResult) apolvalue.getValue();
            np.setApol(apolresult.doubleValue());

            try {
                BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
                bcutDescriptor.initialise(builder);
                np.bcutDescriptor = new ArrayList<>();
                DescriptorValue bcutvalue = bcutDescriptor.calculate(acFull);
                DoubleArrayResult bcutResults = (DoubleArrayResult) bcutvalue.getValue();
                for (int i = 0; i < 6; i++) {
                    np.bcutDescriptor.add(bcutResults.get(i));
                }
            }catch(ArrayIndexOutOfBoundsException e){
                e.printStackTrace();
            }

            BPolDescriptor bPolDescriptor = new BPolDescriptor();
            bPolDescriptor.initialise(builder);
            DescriptorValue bpolvalue = bPolDescriptor.calculate(acFull);
            DoubleResult bpolresult = (DoubleResult) bpolvalue.getValue();
            np.setBpol(bpolresult.doubleValue());

            EccentricConnectivityIndexDescriptor eccentricConnectivityIndexDescriptor = new EccentricConnectivityIndexDescriptor();
            eccentricConnectivityIndexDescriptor.initialise(builder);
            DescriptorValue eccenValue = eccentricConnectivityIndexDescriptor.calculate(acFull);
            IntegerResult eccenResult = (IntegerResult) eccenValue.getValue();
            np.setEccentricConnectivityIndexDescriptor(eccenResult.intValue());

            FMFDescriptor fmfDescriptor = new FMFDescriptor();
            fmfDescriptor.initialise(builder);
            DescriptorValue fmfValue = fmfDescriptor.calculate(acFull);
            DoubleResult fmfResult = (DoubleResult) fmfValue.getValue();
            np.setFmfDescriptor(fmfResult.doubleValue());

            FractionalCSP3Descriptor fractionalCSP3Descriptor = new FractionalCSP3Descriptor();
            fractionalCSP3Descriptor.initialise(builder);
            DescriptorValue fsp3value = fractionalCSP3Descriptor.calculate(acFull);
            DoubleResult fsp3result = (DoubleResult) fsp3value.getValue();
            np.setFsp3(fsp3result.doubleValue());

            FragmentComplexityDescriptor fragmentComplexityDescriptor = new FragmentComplexityDescriptor();
            fragmentComplexityDescriptor.initialise(builder);
            DescriptorValue fcdValue = fragmentComplexityDescriptor.calculate(acFull);
            DoubleResult fcdResults = (DoubleResult) fcdValue.getValue();
            np.setFragmentComplexityDescriptor(fcdResults.doubleValue());


            GravitationalIndexDescriptor gravitationalIndexDescriptor = new GravitationalIndexDescriptor();
            gravitationalIndexDescriptor.initialise(builder);
            DescriptorValue gravValue = gravitationalIndexDescriptor.calculate(acFull);
            DoubleArrayResult gravResults = (DoubleArrayResult) gravValue.getValue();
            np.setGravitationalIndexHeavyAtoms(gravResults.get(0));

            HBondAcceptorCountDescriptor hBondAcceptorCountDescriptor = new HBondAcceptorCountDescriptor();
            hBondAcceptorCountDescriptor.initialise(builder);
            DescriptorValue nhbaccValue = hBondAcceptorCountDescriptor.calculate(acFull);
            IntegerResult nhbaccResult = (IntegerResult) nhbaccValue.getValue();
            np.sethBondAcceptorCount(nhbaccResult.intValue());

            HBondDonorCountDescriptor hBondDonorCountDescriptor = new HBondDonorCountDescriptor();
            hBondDonorCountDescriptor.initialise(builder);
            DescriptorValue nhbdonValue = hBondDonorCountDescriptor.calculate(acFull);
            IntegerResult nhbdonResult = (IntegerResult) nhbdonValue.getValue();
            np.sethBondDonorCount(nhbdonResult.intValue());

            HybridizationRatioDescriptor hybridizationRatioDescriptor = new HybridizationRatioDescriptor();
            hybridizationRatioDescriptor.initialise(builder);
            DescriptorValue hybridRatioValue = hybridizationRatioDescriptor.calculate(acFull);
            DoubleResult hybridRationResult = (DoubleResult) hybridRatioValue.getValue();
            np.setHybridizationRatioDescriptor(hybridRationResult.doubleValue());

            KappaShapeIndicesDescriptor kappaShapeIndicesDescriptor = new KappaShapeIndicesDescriptor();
            kappaShapeIndicesDescriptor.initialise(builder);
            DescriptorValue kappaShapeValues = kappaShapeIndicesDescriptor.calculate(acFull);
            DoubleArrayResult kappaShapeResults = (DoubleArrayResult) kappaShapeValues.getValue();
            np.setKappaShapeIndex1(kappaShapeResults.get(0));
            np.setKappaShapeIndex2(kappaShapeResults.get(1));
            np.setKappaShapeIndex3(kappaShapeResults.get(2));

            MannholdLogPDescriptor mannholdLogPDescriptor = new MannholdLogPDescriptor();
            mannholdLogPDescriptor.initialise(builder);
            DescriptorValue manholdLogpValues = mannholdLogPDescriptor.calculate(acFull);
            DoubleResult manholdLogpResult = (DoubleResult) manholdLogpValues.getValue();
            np.setManholdlogp(manholdLogpResult.doubleValue());

            PetitjeanNumberDescriptor petitjeanNumberDescriptor = new PetitjeanNumberDescriptor();
            petitjeanNumberDescriptor.initialise(builder);
            DescriptorValue petitjeanNumnerValue = petitjeanNumberDescriptor.calculate(acFull);
            DoubleResult petitjeanResult = (DoubleResult) petitjeanNumnerValue.getValue();
            np.setPetitjeanNumber(petitjeanResult.doubleValue());

            PetitjeanShapeIndexDescriptor petitjeanShapeIndexDescriptor = new PetitjeanShapeIndexDescriptor();
            petitjeanShapeIndexDescriptor.initialise(builder);
            DescriptorValue petitjeanShapeValues = petitjeanShapeIndexDescriptor.calculate(acFull);
            DoubleArrayResult petitjeanShapeResults = (DoubleArrayResult) petitjeanShapeValues.getValue();
            np.setPetitjeanShapeTopo(petitjeanShapeResults.get(0));
            np.setPetitjeanShapeGeom(petitjeanShapeResults.get(1));

            RuleOfFiveDescriptor ruleOfFiveDescriptor = new RuleOfFiveDescriptor();
            ruleOfFiveDescriptor.initialise(builder);
            DescriptorValue ruleOfFiveValue = ruleOfFiveDescriptor.calculate(acFull);
            IntegerResult ruleOfFiveResult = (IntegerResult) ruleOfFiveValue.getValue();
            np.setLipinskiRuleOf5Failures(ruleOfFiveResult.intValue());

                /*SmallRingDescriptor smallRingDescriptor = new SmallRingDescriptor();
                smallRingDescriptor.initialise(builder);
                DescriptorValue smallRingValue = smallRingDescriptor.calculate(acFull);
                IntegerResult smallRingResult = (IntegerResult) smallRingValue.getValue();
                np.setNumberSmallRingsDescriptor(smallRingResult.intValue());*/


            SpiroAtomCountDescriptor spiroAtomCountDescriptor = new SpiroAtomCountDescriptor();
            spiroAtomCountDescriptor.initialise(builder);
            DescriptorValue spiroatomValue = spiroAtomCountDescriptor.calculate(acFull);
            IntegerResult spiroatomResult = (IntegerResult) spiroatomValue.getValue();
            np.setNumberSpiroAtoms(spiroatomResult.intValue());

            VABCDescriptor vabcDescriptor = new VABCDescriptor();
            vabcDescriptor.initialise(builder);
            DescriptorValue vabcValue = vabcDescriptor.calculate(acFull);
            DoubleResult vabcResult = (DoubleResult) vabcValue.getValue();
            np.setVabcDescriptor(vabcResult.doubleValue());

            VAdjMaDescriptor vAdjMaDescriptor = new VAdjMaDescriptor();
            vAdjMaDescriptor.initialise(builder);
            DescriptorValue vadjmaValue = vAdjMaDescriptor.calculate(acFull);
            DoubleResult vadjmaResult = (DoubleResult) vadjmaValue.getValue();
            np.setVertexAdjMagnitude(vadjmaResult.doubleValue());

            WienerNumbersDescriptor wienerNumbersDescriptor = new WienerNumbersDescriptor();
            wienerNumbersDescriptor.initialise(builder);
            DescriptorValue wienerNumbersValue = wienerNumbersDescriptor.calculate(acFull);
            DoubleArrayResult wienerNumbersResult = (DoubleArrayResult) wienerNumbersValue.getValue();
            np.setWeinerPathNumber(wienerNumbersResult.get(0));
            np.setWeinerPolarityNumber(wienerNumbersResult.get(1));

            XLogPDescriptor xLogPDescriptor = new XLogPDescriptor();
            xLogPDescriptor.initialise(builder);
            DescriptorValue xlogpValues = xLogPDescriptor.calculate(acFull);
            DoubleResult xlogpResult = (DoubleResult) xlogpValues.getValue();
            np.setXlogp(xlogpResult.doubleValue());

            ZagrebIndexDescriptor zagrebIndexDescriptor = new ZagrebIndexDescriptor();
            zagrebIndexDescriptor.initialise(builder);
            DescriptorValue zagrebIndexValue = zagrebIndexDescriptor.calculate(acFull);
            DoubleResult zagrebIndexresult = (DoubleResult) zagrebIndexValue.getValue();
            np.setZagrebIndex(zagrebIndexresult.doubleValue());

            TPSADescriptor tpsaDescriptor = new TPSADescriptor();
            tpsaDescriptor.initialise(builder);
            DescriptorValue tpsaValue = tpsaDescriptor.calculate(acFull);
            DoubleResult tpsaResult = (DoubleResult) tpsaValue.getValue();
            np.setTopoPSA(tpsaResult.doubleValue());

            FractionalPSADescriptor fractionalPSADescriptor = new FractionalPSADescriptor();
            fractionalPSADescriptor.initialise(builder);
            DescriptorValue ftpsaValue = fractionalPSADescriptor.calculate(acFull);
            DoubleResult ftpsaResult = (DoubleResult) ftpsaValue.getValue();
            np.setTpsaEfficiency(ftpsaResult.doubleValue());



            //compute the clean smiles
            SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique); //Unique - canonical SMILES string, different atom ordering produces the same* SMILES. No isotope or stereochemistry encoded.
            try {
                IAtomContainer nm = AtomContainerManipulator.removeHydrogens(acFull);
                String cleanSmiles =  smilesGenerator.create(nm);
                np.setClean_smiles(cleanSmiles);
            } catch (CDKException e) {
                e.printStackTrace();
            }



            uniqueNaturalProductRepository.save(np);



        } catch (CDKException | OutOfMemoryError e) {
            e.printStackTrace();
        }
        return np;


    }





}
