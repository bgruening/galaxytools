package ${package}.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.regex.Pattern;

import org.snakeyaml.engine.v2.api.Load;
import org.snakeyaml.engine.v2.api.LoadSettings;
import org.snakeyaml.engine.v2.nodes.Tag;
import org.snakeyaml.engine.v2.resolver.ScalarResolver;
import org.snakeyaml.engine.v2.schema.CoreSchema;

public class YamlUtils {

	public static Map<String, Object> mapFromString(final String text) {
		LoadSettings settings = LoadSettings.builder().setSchema(new CoreSchema()).build();
		Load load = new Load(settings);
		final Map<String, Object> result = (Map<String, Object>) load.loadFromString(text);
		return result;
	}

	public static List<Object> listFromString(final String text) {
		LoadSettings settings = LoadSettings.builder().setSchema(new CoreSchema()).build();
		Load load = new Load(settings);
		final List<Object> result = (List<Object>) load.loadFromString(text);
		return result;
	}
}
