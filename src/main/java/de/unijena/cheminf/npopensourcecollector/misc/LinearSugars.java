/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.unijena.cheminf.npopensourcecollector.misc;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
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
 * @author kalai
 */
public class LinearSugars {

    UniversalIsomorphismTester isomorphismTester = null;

    public boolean hasSugarChains(IAtomContainer target, int ringCount) throws CDKException {
        if (ringCount == 0) {
            for (IAtomContainer query : sugarChains) {
                if (isomorphismTester.isSubgraph(target, query)) {
                    return true;
                }
            }
        }
        return false;
    }

    public static class LinearSugarsGeneratorHolder {
        private static final LinearSugars INSTANCE = new LinearSugars();
    }

    public static LinearSugars getInstance() {
        return LinearSugarsGeneratorHolder.INSTANCE;
    }

    private LinearSugars() {
        isomorphismTester = new UniversalIsomorphismTester();
        sugarChains = getSugarChains();
    }

    public static List<IAtomContainer> sugarChains;

    private List<IAtomContainer> getSugarChains() {

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
}
