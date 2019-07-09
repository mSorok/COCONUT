package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.*;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;



@Service
public class MolecularFeaturesComputationService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;

    public void doWork(){
        System.out.println("Calculating additional features for unique molecules");

        //Constructors for descriptors
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();

        List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();

        for(UniqueNaturalProduct np : allNP){
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


                uniqueNaturalProductRepository.save(np);



            } catch (CDKException e) {
                e.printStackTrace();
            }



            uniqueNaturalProductRepository.save(np);
        }


    }
}
