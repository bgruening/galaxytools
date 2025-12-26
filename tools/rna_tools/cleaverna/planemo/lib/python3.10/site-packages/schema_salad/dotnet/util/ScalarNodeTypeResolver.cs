#pragma warning disable CS8769
using System.Text.RegularExpressions;
using YamlDotNet.Core.Events;
using YamlDotNet.Serialization;

namespace ${project_name};
internal class ScalarNodeTypeResolver : INodeTypeResolver
{
    bool INodeTypeResolver.Resolve(NodeEvent nodeEvent, ref Type currentType)
    {
        if (currentType == typeof(object))
        {
            if (nodeEvent is Scalar scalar && !scalar.IsQuotedImplicit)
            {
                // Expressions taken from https://github.com/aaubry/YamlDotNet/blob/feat-schemas/YamlDotNet/Core/Schemas/JsonSchema.cs

                if (Regex.IsMatch(scalar.Value, @"^(true|false|True|False)$", RegexOptions.IgnorePatternWhitespace))
                {
                    currentType = typeof(bool);
                    return true;
                }

                if (Regex.IsMatch(scalar.Value, @"^-? ( 0 | [1-9] [0-9]* )$", RegexOptions.IgnorePatternWhitespace))
                {
                    currentType = typeof(int);
                    return true;
                }

                if (Regex.IsMatch(
                    scalar.Value,
                    @"^-? ( 0 | [1-9] [0-9]* ) ( \. [0-9]* )? ( [eE] [-+]? [0-9]+ )?$",
                    RegexOptions.IgnorePatternWhitespace))
                {
                    currentType = typeof(double);
                    return true;
                }

                // Add more cases here if needed
            }
        }

        return false;
    }
}
