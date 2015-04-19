/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

package bacter.model.pop;

import bacter.ConversionGraph;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Unit tests for BSP implementation.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SkylinePopulationFunctionTest {

    @Test
    public void test() throws Exception {
        String acgString = "[&15,0,1.3905355989030808,31,770,1.597708055397074] " +
                        "[&30,931,2.4351280458424904,36,2486,3.78055549386568] " +
                        "[&15,941,2.0439957300083322,38,2364,6.911056700367016] " +
                        "[&36,1091,4.285505683622974,38,2589,9.867725913197855] " +
                        "((((10:0.5385300170206817,(17:0.116794353049212," +
                        "((3:0.039229346597297564,12:0.039229346597297564)23:0.04582913870888949," +
                        "13:0.08505848530618705)24:0.03173586774302495)26:0.4217356639714697)28:1.8114199763246093," +
                        "((8:0.10883006062265468,2:0.10883006062265468)25:0.556428062025291," +
                        "(6:0.5393311342677402,11:0.5393311342677402)29:0.12592698838020555)31:1.6846918706973453)34:1.4536824928125807," +
                        "(1:0.47184545557390367,14:0.47184545557390367)27:3.331787030583968)37:2.9704369411362554," +
                        "(((15:2.0624287390593707,((16:0.01825347077733299,19:0.01825347077733299)21:0.7668749128372041," +
                        "(7:0.008018731329538273,9:0.008018731329538273)20:0.7771096522849988)32:1.2773003554448337)33:0.7487092404613747," +
                        "4:2.8111379795207454)35:0.1331794525400949,((0:0.0243537216663141," +
                        "5:0.0243537216663141)22:0.5681537100482162,18:0.5925074317145304)30:2.35181000034631)36:3.829751995233287)38:0.0";

        ConversionGraph acg = new ConversionGraph();
        acg.initByName("sequenceLength", 10000, "fromString", acgString);

        SkylinePopulationFunction skyline = new SkylinePopulationFunction();
        skyline.initByName(
                "acg", acg,
                "popSizes", new RealParameter("5.0 1.0 5.0 1.0"),
                "groupSizes", new IntegerParameter("0 0 0 0"));

        for (double t = 0.0; t<10; t += 0.01)
            assertTrue(Math.abs(t-skyline.getInverseIntensity(skyline.getIntensity(t)))<1e-14);
    }
}
