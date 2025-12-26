using ${project_name};
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test;

[TestClass]
public class UtilitiesTests
{
    [TestMethod]
    public void TestShortname()
    {
        Assert.AreEqual("homo_sapiens", Utilities.Shortname("file:///Users/jdidion/projects/cwlScala/target/test-classes/CommandLineTools/conformance/#anon_enum_inside_array_inside_schemadef.cwl/first/user_type_2/species/homo_sapiens"));
        Assert.AreEqual("GRCh37", Utilities.Shortname("file:///home/michael/cwljava/src/test/resources/org/w3id/cwl/cwl1_2/utils/valid_anon_enum_inside_array_inside_schemadef.cwl#vcf2maf_params/ncbi_build/GRCh37"));
        Assert.AreEqual("foo", Utilities.Shortname("http://example.com/foo"));
        Assert.AreEqual("bar", Utilities.Shortname("http://example.com/#bar"));
        Assert.AreEqual("bar", Utilities.Shortname("http://example.com/foo/bar"));
        Assert.AreEqual("bar", Utilities.Shortname("http://example.com/foo#bar"));
        Assert.AreEqual("bar", Utilities.Shortname("http://example.com/#foo/bar"));
        Assert.AreEqual("baz", Utilities.Shortname("http://example.com/foo#bar/baz"));
    }
}
