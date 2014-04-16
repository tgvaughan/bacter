/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package argbeast;

import java.util.List;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * Unit test for marginal tree traversal.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalTreeTest {
    
    public MarginalTreeTest() {
    }

    @Test
    public void test() throws Exception {

        // Recombination graph
        String str = "[&2,2759,0.3260126313706676,10,2808,0.42839862922656696] "
                + "[&10,6692,0.3381366423491633,2,6693,0.5683827224649434] "
                + "[&10,8069,0.2807615297583804,14,8160,0.3415740002783274] "
                + "(((0:0.04916909893812008,1:0.04916909893812008)10:0.5465237639426681,"
                + "(4:0.3773111326866937,(((8:0.22180790639747835,"
                + "(3:0.07561592852503513,6:0.07561592852503513)11:0.14619197787244323)"
                + "13:0.010206467073885589,9:0.23201437347136394)14:0.116542689187905,"
                + "(7:0.10746702934931932,5:0.10746702934931932)12:0.24109003330994963)"
                + "15:0.02875407002742475)16:0.21838173019409446)17:1.1073878800617445,"
                + "2:1.7030807429425328)18:0.0";
        
        RecombinationGraph arg = new RecombinationGraph();
        arg.initByName("fromString", str, "sequenceLength", 10000);
        List<Recombination> recombs = arg.getRecombinations();
        
        for (Recombination recomb : arg.getRecombinations())
            System.out.println(arg.getMarginalNewick(recomb));
        
        // Test root nodes values
        assertEquals(18, arg.getMarginalRoot(recombs.get(0)).getNr());
        assertEquals(17, arg.getMarginalRoot(recombs.get(1)).getNr());
        assertEquals(18, arg.getMarginalRoot(recombs.get(2)).getNr());
        assertEquals(18, arg.getMarginalRoot(recombs.get(3)).getNr());
        
        // Test root heights
        assertTrue(Math.abs(arg.getMarginalNodeHeight(arg.getMarginalRoot(
                recombs.get(0)), recombs.get(0)) - 1.7031)<1e-3);
        assertTrue(Math.abs(arg.getMarginalNodeHeight(arg.getMarginalRoot(
                recombs.get(1)), recombs.get(1)) - 0.59569)<1e-3);
        assertTrue(Math.abs(arg.getMarginalNodeHeight(arg.getMarginalRoot(
                recombs.get(2)), recombs.get(2)) - 1.7031)<1e-3);
        assertTrue(Math.abs(arg.getMarginalNodeHeight(arg.getMarginalRoot(
                recombs.get(3)), recombs.get(3)) - 1.7031)<1e-3);
    }
}
