package test.commands;

import static org.junit.Assert.assertEquals;
import junit.framework.JUnit4TestAdapter;

import org.junit.Test;

import commands.CommandRouter;

public class TestCommandRouter {
	
	public static final String[] threeItems = {"1","2","3"};
	public static final String[] firstItemGone = {"2","3"};

	@Test public void testCreateCmdArgs() {
		assertEquals(firstItemGone, CommandRouter.createCmdArgs(threeItems));
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestCommandRouter.class);
	}

}
