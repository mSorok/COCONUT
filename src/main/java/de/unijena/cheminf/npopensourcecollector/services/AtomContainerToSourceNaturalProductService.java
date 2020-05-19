package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.springframework.stereotype.Service;

@Service
public class AtomContainerToSourceNaturalProductService {



    public SourceNaturalProduct createSNPlInstance(IAtomContainer ac) {
        SourceNaturalProduct np = new SourceNaturalProduct();

        np.setSource(ac.getProperty("SOURCE"));

        np.setIdInSource(ac.getID());

        np.setOriginalInchi(ac.getProperty("ORIGINAL_INCHI"));
        np.setOriginalInchiKey(ac.getProperty("ORIGINAL_INCHIKEY"));
        np.setOriginalSmiles(ac.getProperty("ORIGINAL_SMILES"));

        np.setSimpleInchi(ac.getProperty("SIMPLE_INCHI"));
        np.setSimpleInchiKey(ac.getProperty("SIMPLE_INCHIKEY"));
        np.setSimpleSmiles(ac.getProperty("SIMPLE_SMILES"));

        if(ac.getProperties().containsKey("ABSOLUTE_SMILES")) {
            np.setAbsoluteSmiles(ac.getProperty("ABSOLUTE_SMILES"));
        }

        if(ac.getProperties().containsKey("CAS")) {
            np.setCas(ac.getProperty("CAS"));
        }


        np.setTotalAtomNumber(ac.getAtomCount());

        int heavyAtomCount = 0;
        for(IAtom a : ac.atoms()){
            if(!a.getSymbol().equals("H")){
                heavyAtomCount=heavyAtomCount+1;
            }
        }
        np.setHeavyAtomNumber(heavyAtomCount);

        np.setAcquisition_date(ac.getProperty("ACQUISITION_DATE"));

        return np;
    }


    public IAtomContainer createAtomContainer(SourceNaturalProduct snp){
        IAtomContainer ac = null;

        try {
            SmilesParser sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
            ac   = sp.parseSmiles( snp.getSimpleSmiles() );
        } catch (InvalidSmilesException e) {
            System.err.println(e.getMessage());
        }

        return ac;
    }

}
