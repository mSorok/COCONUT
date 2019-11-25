/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.unijena.cheminf.npopensourcecollector.misc;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.List;

/**
 * Singleton class to represent linear sugar chains that are to be removed from input molecule
 * during sugar curation. The set here is limited to the author's knowledge and more chains can be
 * added.
 *
 * @author kalai, msorok
 */
public class LinearSugars {

    UniversalIsomorphismTester isomorphismTester = null;

    public static List<IAtomContainer> sugarChains;

    public List<DfPattern> patternList;

    public IAtomContainer removeLinearSugars(IAtomContainer target) throws CDKException{
        IAtomContainer newMolecule = null;
        try {
            newMolecule = target.clone();

            for(IAtomContainer sugarChain : sugarChains){
                if( isomorphismTester.isSubgraph(target, sugarChain)){

                    //remove sugar chain
                }
            }

        } catch (CloneNotSupportedException e) {
            return null;
        }


        return newMolecule;

    }



    public static class LinearSugarsGeneratorHolder {
        private static final LinearSugars INSTANCE = new LinearSugars();
    }

    public static LinearSugars getInstance() {
        return LinearSugarsGeneratorHolder.INSTANCE;
    }

    public LinearSugars() {
        sugarChains = getSugarChains();
        patternList = getSugarPatterns();
    }


    public List<IAtomContainer> getSugarChains() {

        List<IAtomContainer> linearSugarChains = new ArrayList<IAtomContainer>();

        String[] smilesList = {"C(C(C(C(C(C=O)O)O)O)O)O", "C(C(CC(C(CO)O)O)O)(O)=O", "C(C(C(CC(=O)O)O)O)O",
                "C(C(C(C(C(CO)O)O)O)=O)O", "C(C(C(C(C(CO)O)O)O)O)O", "C(C(C(C(CC=O)O)O)O)O", "occ(o)co",
                "OCC(O)C(O)C(O)C(O)CO", "O=CC(O)C(O)C(O)C(O)CO", "CC(=O)OCC(O)CO", "CCCCC(O)C(=O)O", "CC(=O)CC(=O)CCC(=O)O",
                "CC(O)C(O)C(=O)O", "O=C(O)CC(O)CC(=O)O", "O=C(O)C(=O)C(=O)C(O)C(O)CO", "CC(O)CC(=O)O", "CC(CCC(=O)O)CC(=O)O",
                "O=C(O)CCC(O)C(=O)O", "O=CC(O)C(O)C(O)C(O)CO", "O=C(CO)C(O)C(O)CO"};
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

        for (String smiles : smilesList) {
            try {
                linearSugarChains.add(sp.parseSmiles(smiles));
            } catch (InvalidSmilesException ex) {
                System.out.println("could not remove linear sugar");
            }
        }
        return linearSugarChains;
    }


    public List<DfPattern> getSugarPatterns(){
        List<DfPattern> patternList = new ArrayList<>();
        for(IAtomContainer sugarAC : sugarChains){
            patternList.add(DfPattern.findSubstructure(sugarAC));

        }

        return patternList;
    }
}
