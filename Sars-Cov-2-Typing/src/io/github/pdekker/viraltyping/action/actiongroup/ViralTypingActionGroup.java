package io.github.pdekker.viraltyping.action.actiongroup;

import com.clcbio.api.free.actions.framework.ActionGroup;
import com.clcbio.api.free.actions.framework.ActionGroup.GroupType;
import com.clcbio.api.free.actions.framework.ActionGroupPlugin;
import com.clcbio.api.free.actions.framework.StaticActionGroupDefinitions;
import com.clcbio.api.free.gui.icon.ClcIcon;
import com.clcbio.api.free.gui.icon.DefaultClcIcon;

public class ViralTypingActionGroup extends ActionGroupPlugin {
	public static final String PLUGIN_GROUP = "free";
	public static final String CLASS_KEY = "ViralTypingActionGroup";

	@Override
	public String getName() {
		return "SARS-CoV-2 Analysis";
	}

	@Override
	public ActionGroup getParent() {
		return StaticActionGroupDefinitions.TOOLBOX_TOP_GROUP;
	}

	@Override
	public ClcIcon getOpenIcon() {
		return new DefaultClcIcon("icons/virus");
	}

	@Override
	public ClcIcon getCloseIcon() {
		return new DefaultClcIcon("icons/virus");
	}

	@Override
	public GroupType getType() {
		return ActionGroup.SUBMENUTYPE;
	}

	@Override
	public String getClassKey() {
		return CLASS_KEY;
	}

	@Override
	public double getVersion() {
		return 0.1;
	}

	@Override
	public int getPreferredMenuLocation() {
		return 0;
	}

}
