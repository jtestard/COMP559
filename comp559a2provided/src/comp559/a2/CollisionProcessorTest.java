package comp559.a2;

import static org.junit.Assert.*;
import junit.framework.Assert;

import org.junit.Test;

public class CollisionProcessorTest {

	@Test
	public void test() {
		CollisionProcessor pr = new CollisionProcessor(null);
		pr.distances[0] = 2;
		pr.distances[1] = 3;
		pr.distances[2] = 4;
		pr.distances[3] = 1;
		Assert.assertEquals(0, pr.smallestDistance(true));
		Assert.assertEquals(3, pr.smallestDistance(false));
	}

}
